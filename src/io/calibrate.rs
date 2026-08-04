//! Choosing a gzip decoder from measured input characteristics.
//!
//! Reusable across tools that map FASTX through this crate — piscem-rs itself,
//! and downstream consumers such as salmon, which query through
//! [`crate::mapping::streaming_query`] and face the same decision. Nothing here
//! depends on the mapping kernel: the probe takes a closure, so the caller
//! supplies whatever "map one record" means for it.
//!
//! # The decision
//!
//! The serial path (niffler/libz-rs) opens exactly one inflate stream per input
//! file. Its supply is `files x per-stream rate` and does **not** respond to
//! `-t`. Demand is `map_threads x per-thread consumption`. The parallel decoder
//! is worth its overhead when demand exceeds supply.
//!
//! Per-thread consumption is the hard term: measured 0.064 GB/s on a 96.3 M
//! k-mer transcriptome versus 0.43 GB/s on a 1.5 M k-mer probe panel, a 6.7x
//! spread. No constant spans that, and two data points do not justify fitting a
//! curve against index size, so it is *measured* rather than predicted.
//!
//! # Two tiers
//!
//! 1. [`forced_choice`] — decisions that follow from the input alone, with no
//!    measurement. These are logical forcings, not heuristics: cases where the
//!    answer is the same whatever the mapping rate turns out to be.
//! 2. [`probe`] + [`decide`] — when nothing forces the answer, measure both
//!    rates from a small prefix and compare. Measured cost of the probe: 0.04 s
//!    for 5,000 reads against a 96.3 M k-mer index, 0.02 s against a probe
//!    panel, single-threaded.
//!
//! # Why probing does not deadlock on a stream
//!
//! Probing means reading a prefix and then re-reading from the start, which a
//! pipe cannot do. That is not a limitation here, because a non-seekable input
//! is already a forced case: `rapidgzip` requires positional reads and falls
//! back to *sequential* decoding when the input is not a regular file, so the
//! parallel path has nothing to offer and the probe is never reached.

use std::path::Path;

/// Decompressed GB/s obtainable from one serial inflate stream.
///
/// Deliberately the conservative (multi-member) figure: measured 2.40 GB/s on
/// single-member archives and 1.45 GB/s on multi-member ones. Underestimating
/// supply biases toward the parallel decoder, which is the cheap error — see
/// [`decide`].
const SERIAL_GBPS_PER_STREAM: f64 = 1.45;

/// What kind of input this is, which decides whether a restart is even possible.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InputKind {
    /// A regular file: supports positional reads, and can be re-opened.
    Regular,
    /// A pipe, FIFO, socket or character device. Cannot be re-read, and
    /// `rapidgzip` degrades to sequential decoding on it.
    Stream,
}

/// Classify `path`. Anything that is not a regular file is a [`InputKind::Stream`].
///
/// Mirrors `rapidgzip`'s own `supports_positional_reads`, which gates its
/// parallel path on `file_type().is_file()`. Keeping the same test here means
/// this crate never hands the parallel decoder an input it would silently
/// decode sequentially.
pub fn classify_input(path: &Path) -> InputKind {
    match std::fs::metadata(path) {
        Ok(md) if md.file_type().is_file() => InputKind::Regular,
        // Unreadable metadata is treated as a stream: the conservative answer,
        // since the serial path works on everything.
        _ => InputKind::Stream,
    }
}

/// Which decoder to use, and why.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Decision {
    pub parallel: bool,
    pub reason: Reason,
}

/// Why a decision was reached. Carried so callers can log it; users hitting an
/// unexpectedly slow run should be able to see which branch fired.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Reason {
    /// Input is not a regular file. The parallel decoder cannot use positional
    /// reads and would decode sequentially anyway, while costing coordinator
    /// threads. Also the case where a probe could not restart the input.
    NonSeekableInput,
    /// Only one mapping thread, so decode can never be the constraint.
    SingleThreaded,
    /// Serial streams already outnumber what the mapping threads can consume by
    /// so much that no mapping rate makes the parallel decoder worthwhile.
    AmpleSerialSupply,
    /// Measured demand exceeds measured supply.
    MeasuredDecodeBound,
    /// Measured supply covers measured demand.
    MeasuredConsumerBound,
}

/// Decisions that follow from the input alone, before any measurement.
///
/// Returns `None` when the answer genuinely depends on how fast mapping
/// consumes data, in which case the caller should [`probe`] and [`decide`].
///
/// These are forcings rather than guesses:
///
/// * A non-regular file cannot be decoded in parallel at all.
/// * With one mapping thread there is nothing to starve.
/// * Below [`MIN_THREADS_PER_FILE`] threads per file, even a free decoder
///   cannot help: serial supply already exceeds what that few threads can
///   consume at the *fastest* per-thread rate ever measured here (0.43 GB/s),
///   so no probe can change the answer.
pub fn forced_choice(kind: InputKind, num_files: usize, map_threads: usize) -> Option<Decision> {
    if kind == InputKind::Stream {
        return Some(Decision { parallel: false, reason: Reason::NonSeekableInput });
    }
    if map_threads <= 1 {
        return Some(Decision { parallel: false, reason: Reason::SingleThreaded });
    }
    if map_threads / num_files.max(1) < MIN_THREADS_PER_FILE {
        return Some(Decision { parallel: false, reason: Reason::AmpleSerialSupply });
    }
    None
}

/// Threads per file below which serial supply is sufficient at any mapping rate.
///
/// `SERIAL_GBPS_PER_STREAM / fastest observed per-thread consumption`
/// = 1.45 / 0.43 ~= 3.4, rounded down. Independently corroborated by the
/// measured crossover surface, where 4 threads/file sits at 1.07-1.09x and 2
/// threads/file at 0.92-1.00x.
const MIN_THREADS_PER_FILE: usize = 3;

/// Rates measured from a prefix of the input.
#[derive(Debug, Clone, Copy)]
pub struct InputRates {
    /// Decompressed GB/s from the one serial inflate stream that produced them.
    pub decode_gbps_per_stream: f64,
    /// Decompressed GB/s a single mapping thread consumes.
    pub map_gbps_per_thread: f64,
    pub records_sampled: usize,
    pub bytes_sampled: u64,
}

/// Measure decode and mapping rates from the first `sample` records.
///
/// `next_record` yields decompressed record bytes from a *serial* reader — the
/// one decoder available before any decision has been made. `map_one` performs
/// whatever the caller considers mapping; it is timed but its result is
/// ignored.
///
/// Both timings come from the same prefix, so a slow disk inflates the decode
/// rate's denominator and not the mapping rate's, which is the correct
/// attribution: decode competes with I/O, mapping does not.
///
/// Returns `None` if the input ended before yielding any records, in which case
/// there is nothing to decide.
pub fn probe<N, M>(sample: usize, mut next_record: N, mut map_one: M) -> Option<InputRates>
where
    N: FnMut() -> Option<u64>,
    M: FnMut(),
{
    use std::time::Instant;

    let mut bytes = 0u64;
    let mut records = 0usize;
    let mut map_elapsed = std::time::Duration::ZERO;
    let decode_start = Instant::now();

    while records < sample {
        let Some(len) = next_record() else { break };
        bytes += len;
        records += 1;
        let t = Instant::now();
        map_one();
        map_elapsed += t.elapsed();
    }

    if records == 0 || bytes == 0 {
        return None;
    }

    // Decode time is wall time spent pulling bytes, i.e. total minus mapping.
    let total = decode_start.elapsed();
    let decode_elapsed = total.saturating_sub(map_elapsed);
    let gb = bytes as f64 / 1e9;

    Some(InputRates {
        decode_gbps_per_stream: safe_rate(gb, decode_elapsed.as_secs_f64()),
        map_gbps_per_thread: safe_rate(gb, map_elapsed.as_secs_f64()),
        records_sampled: records,
        bytes_sampled: bytes,
    })
}

/// Rate that degrades to a sentinel rather than dividing by ~zero.
fn safe_rate(gb: f64, secs: f64) -> f64 {
    if secs <= 1e-9 { f64::INFINITY } else { gb / secs }
}

/// Compare measured demand against measured supply.
///
/// `margin` biases the comparison, and should be **below** 1.0. The payoff is
/// sharply asymmetric: enabling the parallel decoder when it is not needed
/// measured -4.9% wall and +2.1% CPU on a transcriptome where decode never
/// binds, while failing to enable it when it is needed costs up to 3x wall. So
/// the tie is broken toward parallel.
pub fn decide(
    rates: &InputRates,
    num_files: usize,
    map_threads: usize,
    margin: f64,
) -> Decision {
    // A probe run on one thread cannot see the memory contention of a full run,
    // so `map_gbps_per_thread` overstates what each thread will really consume,
    // and the supply figure is the conservative multi-member one. Both biases
    // push toward the parallel decoder, which is the cheap error.
    let supply = num_files.max(1) as f64 * rates.decode_gbps_per_stream.min(SERIAL_GBPS_PER_STREAM);
    let demand = map_threads as f64 * rates.map_gbps_per_thread;

    if demand > supply * margin {
        Decision { parallel: true, reason: Reason::MeasuredDecodeBound }
    } else {
        Decision { parallel: false, reason: Reason::MeasuredConsumerBound }
    }
}

/// Default bias for [`decide`].
pub const DEFAULT_MARGIN: f64 = 0.8;

/// Records to sample when probing.
///
/// Measured cost single-threaded: 0.04 s against a 96.3 M k-mer transcriptome,
/// 0.02 s against a 1.5 M k-mer probe panel. Small enough to run
/// unconditionally whenever the answer is not forced.
pub const DEFAULT_PROBE_RECORDS: usize = 5_000;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn streams_never_take_the_parallel_path() {
        // The FIFO case: rapidgzip would decode it sequentially anyway, and a
        // probe could not rewind it.
        let d = forced_choice(InputKind::Stream, 1, 64).unwrap();
        assert!(!d.parallel);
        assert_eq!(d.reason, Reason::NonSeekableInput);
    }

    #[test]
    fn single_thread_is_forced_serial() {
        let d = forced_choice(InputKind::Regular, 1, 1).unwrap();
        assert!(!d.parallel);
        assert_eq!(d.reason, Reason::SingleThreaded);
    }

    #[test]
    fn many_files_and_few_threads_is_forced_serial() {
        // 8 files, 16 threads = 2 threads/file: measured 0.92-1.00x, and no
        // mapping rate can make it pay.
        assert_eq!(
            forced_choice(InputKind::Regular, 8, 16).unwrap().reason,
            Reason::AmpleSerialSupply
        );
    }

    #[test]
    fn ambiguous_cases_defer_to_measurement() {
        // 2 files, 32 threads = 16 threads/file. Whether this pays depends
        // entirely on how fast mapping consumes, so it must not be forced.
        assert!(forced_choice(InputKind::Regular, 2, 32).is_none());
        // 1 file, 32 threads: same.
        assert!(forced_choice(InputKind::Regular, 1, 32).is_none());
    }

    /// The two regimes actually measured, fed through `decide` as rates.
    #[test]
    fn decide_matches_the_measured_regimes() {
        // Probe panel, 2 files, 32 threads: demand 32 x 0.43 = 13.8 GB/s
        // against supply 2 x 1.45 = 2.9. Measured 3.05x for parallel.
        let panel = InputRates {
            decode_gbps_per_stream: 1.45,
            map_gbps_per_thread: 0.43,
            records_sampled: 5000,
            bytes_sampled: 1 << 21,
        };
        assert!(decide(&panel, 2, 32, DEFAULT_MARGIN).parallel);

        // 96.3M-k-mer transcriptome, 1 file, 32 threads: demand 32 x 0.064 =
        // 2.05 against supply 1.45. Close, and the asymmetry breaks it toward
        // parallel -- which measured -4.9% wall / +2.1% CPU, an acceptable loss.
        let txome = InputRates {
            decode_gbps_per_stream: 1.45,
            map_gbps_per_thread: 0.064,
            records_sampled: 5000,
            bytes_sampled: 1 << 21,
        };
        assert!(decide(&txome, 1, 32, DEFAULT_MARGIN).parallel);

        // Same index, but eight files: supply 11.6 swamps demand 2.05.
        assert!(!decide(&txome, 8, 32, DEFAULT_MARGIN).parallel);
    }

    #[test]
    fn probe_reports_none_on_empty_input() {
        assert!(probe(100, || None, || {}).is_none());
    }

    #[test]
    fn probe_attributes_time_to_the_right_rate() {
        let mut left = 10;
        let rates = probe(
            10,
            move || {
                if left == 0 {
                    return None;
                }
                left -= 1;
                std::thread::sleep(std::time::Duration::from_micros(200));
                Some(1_000_000)
            },
            || std::thread::sleep(std::time::Duration::from_micros(800)),
        )
        .unwrap();
        assert_eq!(rates.records_sampled, 10);
        assert_eq!(rates.bytes_sampled, 10_000_000);
        // Mapping got ~4x the time, so its rate must be the lower of the two.
        assert!(
            rates.map_gbps_per_thread < rates.decode_gbps_per_stream,
            "map {} decode {}",
            rates.map_gbps_per_thread,
            rates.decode_gbps_per_stream
        );
    }
}
