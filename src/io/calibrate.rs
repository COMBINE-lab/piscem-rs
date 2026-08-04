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
/// A ceiling on the measured producer rate, not an estimate of it.
///
/// Whole-file serial decode measures 2.24 GB/s (single-member) and 2.49 GB/s
/// (multi-member) on this hardware, cold and warm alike — page-cache state made
/// no measurable difference. The producer also parses, so its real rate is
/// below this; the constant exists only to stop an anomalous probe from
/// inventing supply that no stream can deliver.
///
/// Set to the lower of the two measured ceilings. Underestimating supply biases
/// toward the parallel decoder, which is the cheap error — see [`decide`].
const SERIAL_GBPS_PER_STREAM: f64 = 2.24;

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
    /// The user asked for this decoder explicitly.
    UserRequested,
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
/// Runs in two phases rather than interleaving them: decode the whole sample
/// into memory, timing that, then map all of it, timing that. Alternating the
/// two on one thread — as this used to — never lets decode run ahead the way it
/// does in the real pipeline, so it measures "decode while also mapping" and
/// understated stream throughput by roughly 1.7x even after discarding a
/// warm-up prefix. Separating the phases removes that coupling: each rate is
/// measured with the other idle, which is the quantity the supply/demand
/// comparison actually wants.
///
/// The buffer is bounded by `sample` records — about 1.7 MB at the default,
/// which is why holding the whole thing is preferable to being clever.
///
/// `next_record` yields decompressed record sequences from a *serial* reader,
/// the one decoder available before any decision has been made. Returns `None`
/// if the input yielded no records.
pub fn probe<N, M>(sample: usize, mut next_record: N, mut map_one: M) -> Option<InputRates>
where
    N: FnMut() -> Option<Vec<u8>>,
    M: FnMut(&[u8]),
{
    use std::time::Instant;

    // Warm-up: opening the file, faulting in buffers and decoding the first
    // block say nothing about steady state. Pulled but not timed.
    let warmup = (sample / 8).clamp(1, WARMUP_RECORDS_MAX);
    for _ in 0..warmup {
        if next_record().is_none() {
            return None;
        }
    }

    // Phase 1: decode, with nothing else competing.
    let mut records: Vec<Vec<u8>> = Vec::with_capacity(sample.saturating_sub(warmup));
    let mut bytes = 0u64;
    let decode_start = Instant::now();
    while records.len() + warmup < sample {
        let Some(rec) = next_record() else { break };
        bytes += rec.len() as u64;
        records.push(rec);
    }
    let decode_elapsed = decode_start.elapsed();

    if records.is_empty() || bytes == 0 {
        return None;
    }

    // Phase 2: map, with decode finished.
    let map_start = Instant::now();
    for rec in &records {
        map_one(rec);
    }
    let map_elapsed = map_start.elapsed();

    let gb = bytes as f64 / 1e9;
    Some(InputRates {
        decode_gbps_per_stream: safe_rate(gb, decode_elapsed.as_secs_f64()),
        map_gbps_per_thread: safe_rate(gb, map_elapsed.as_secs_f64()),
        records_sampled: records.len(),
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

/// A user's explicit instruction about decoder selection.
///
/// The escape hatch for someone who knows their machine better than a 5,000
/// record probe does — an NFS mount where decode is unusually slow, a shared
/// node where spending cores on decode is antisocial, or simply reproducing a
/// measurement. `Auto` is the default and nothing changes for users who do not
/// reach for it.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum DecoderPreference {
    /// Decide by forcing rules, then measurement.
    #[default]
    Auto,
    /// Always use the serial (libz-rs) path.
    Serial,
    /// Use the parallel decoder wherever it is possible at all, with an
    /// optional per-file worker ceiling.
    ///
    /// Still yields to [`Reason::NonSeekableInput`]: on a pipe the parallel
    /// decoder degrades to sequential decoding, so honouring the request
    /// literally would spend coordinator threads to obtain nothing. A
    /// preference cannot make an input seekable.
    Parallel { workers_per_file: Option<usize> },
}

impl DecoderPreference {
    /// Resolve `s` as given on a command line: `auto`, `serial`, `parallel`, or
    /// `parallel=N`.
    pub fn parse(s: &str) -> Result<Self, String> {
        let s = s.trim();
        let (head, tail) = match s.split_once('=') {
            Some((h, t)) => (h, Some(t)),
            None => (s, None),
        };
        match (head.to_ascii_lowercase().as_str(), tail) {
            ("auto", None) => Ok(Self::Auto),
            ("serial" | "libz" | "libdeflate", None) => Ok(Self::Serial),
            ("parallel" | "rapidgzip", None) => Ok(Self::Parallel { workers_per_file: None }),
            ("parallel" | "rapidgzip", Some(n)) => n
                .parse::<usize>()
                .map(|w| Self::Parallel { workers_per_file: Some(w.max(1)) })
                .map_err(|_| format!("`{n}` is not a worker count")),
            _ => Err(format!(
                "expected `auto`, `serial`, `parallel`, or `parallel=N`, got `{s}`"
            )),
        }
    }
}

/// How the input is compressed, as far as the decoder choice cares.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InputCompression {
    /// gzip — the only format either decoder path competes over.
    Gzip,
    /// Plain text, or a codec (zstd, bz2) that only the serial path handles.
    NotGzip,
    /// Not determined, because determining it would consume the input.
    Unknown,
}

/// Sniff whether `path` is gzip.
///
/// Only ever called for [`InputKind::Regular`]. Sniffing a stream would consume
/// the bytes it read — exactly the double-open that made FIFO input hang — so a
/// non-regular input reports [`InputCompression::Unknown`] rather than being
/// opened speculatively.
pub fn detect_compression(path: &Path, kind: InputKind) -> InputCompression {
    use std::io::Read;
    if kind != InputKind::Regular {
        return InputCompression::Unknown;
    }
    let Ok(mut f) = std::fs::File::open(path) else {
        return InputCompression::Unknown;
    };
    let mut magic = [0u8; 2];
    match f.read_exact(&mut magic) {
        Ok(()) if magic == [0x1f, 0x8b] => InputCompression::Gzip,
        Ok(()) => InputCompression::NotGzip,
        Err(_) => InputCompression::Unknown,
    }
}

/// A request that cannot be carried out as asked.
///
/// Surfaced so the user is told their flag was overridden, rather than quietly
/// getting something else. Someone who passes `--decoder parallel` and sees no
/// speedup deserves to know why.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PreferenceConflict {
    /// `parallel` asked for on an input that cannot support positional reads.
    ParallelOnStream,
    /// A decoder named for an input neither decoder path applies to.
    NotGzipInput,
}

impl std::fmt::Display for PreferenceConflict {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::ParallelOnStream => write!(
                f,
                "--decoder parallel was requested, but the input is not a regular file \
                 (a pipe, FIFO or process substitution). The parallel decoder needs \
                 positional reads and would decode sequentially anyway, so the serial \
                 decoder is being used instead. To get parallel decoding, pass a real \
                 file rather than a stream."
            ),
            Self::NotGzipInput => write!(
                f,
                "--decoder was set, but the input is not gzip-compressed, so the setting \
                 has no effect: only gzip input has both a serial and a parallel decoder. \
                 Plain, zstd and bzip2 inputs are always read by the serial path."
            ),
        }
    }
}

/// Report any way in which `pref` cannot be honoured for this input.
pub fn preference_conflict(
    pref: DecoderPreference,
    kind: InputKind,
    compression: InputCompression,
) -> Option<PreferenceConflict> {
    if pref == DecoderPreference::Auto {
        return None;
    }
    if compression == InputCompression::NotGzip {
        return Some(PreferenceConflict::NotGzipInput);
    }
    if matches!(pref, DecoderPreference::Parallel { .. }) && kind == InputKind::Stream {
        return Some(PreferenceConflict::ParallelOnStream);
    }
    None
}

/// Apply an explicit preference, if it settles the question.
///
/// Returns `None` for [`DecoderPreference::Auto`], and for a parallel request
/// on an input that cannot support it, so the caller falls through to the
/// normal path in both cases.
pub fn preference_choice(pref: DecoderPreference, kind: InputKind) -> Option<Decision> {
    match pref {
        DecoderPreference::Auto => None,
        DecoderPreference::Serial => {
            Some(Decision { parallel: false, reason: Reason::UserRequested })
        }
        DecoderPreference::Parallel { .. } if kind == InputKind::Stream => {
            Some(Decision { parallel: false, reason: Reason::NonSeekableInput })
        }
        DecoderPreference::Parallel { .. } => {
            Some(Decision { parallel: true, reason: Reason::UserRequested })
        }
    }
}

/// Default bias for [`decide`].
pub const DEFAULT_MARGIN: f64 = 0.8;

/// Upper bound on the untimed warm-up prefix. See [`probe`].
const WARMUP_RECORDS_MAX: usize = 1_000;

/// Records to sample when probing.
///
/// Measured cost single-threaded: 0.04 s against a 96.3 M k-mer transcriptome,
/// 0.02 s against a 1.5 M k-mer probe panel. Small enough to run
/// unconditionally whenever the answer is not forced.
pub const DEFAULT_PROBE_RECORDS: usize = 5_000;

/// Measure the producer and consumer sides of the pipeline on a prefix.
///
/// The split follows `paraseq`'s own: `RecordSet::fill` performs decode **and**
/// parsing under the per-file mutex, so it is the *producer* and is serialized
/// per input file; mapping the filled batch is the *consumer* and runs on every
/// thread. Supply is therefore decode+parse throughput, not raw inflate.
///
/// That distinction matters and is easy to get backwards. Raw serial inflate
/// measures 2.2-2.5 GB/s on these inputs, but a probe reporting that as supply
/// would overstate it by whatever parsing costs, biasing toward the serial
/// decoder -- the expensive error. Parsing has to be charged to the producer
/// because that is where the real pipeline pays it.
///
/// Record sets are kept whole and mapped in place, matching the zero-copy path
/// the processors take. An earlier version copied each record out with
/// `to_vec`, charging thousands of allocations to the producer that the real
/// pipeline never performs.
pub fn probe_path<M>(
    path: &std::path::Path,
    sample: usize,
    mut map_one: M,
) -> anyhow::Result<Option<InputRates>>
where
    M: FnMut(&[u8]),
{
    use paraseq::fastx;
    use paraseq::prelude::*;
    use std::time::Instant;

    // Serial open: no decoder threads, works on every input kind.
    let (reader, _fmt) = niffler::send::from_path(path)?;
    let mut reader = fastx::Reader::new(reader)?;

    // Warm-up: one set, pulled and discarded, so opening the file and decoding
    // the first block are not charged to steady state.
    let mut warm = reader.new_record_set();
    if !warm.fill(&mut reader).unwrap_or(false) {
        return Ok(None);
    }
    drop(warm);

    // Producer: decode + parse, nothing else running.
    let mut sets: Vec<fastx::RecordSet> = Vec::new();
    let mut records = 0usize;
    let mut bytes = 0u64;
    let produce_start = Instant::now();
    while records < sample {
        let mut rs = reader.new_record_set();
        match rs.fill(&mut reader) {
            Ok(true) => {}
            _ => break,
        }
        for rec in rs.iter().filter_map(Result::ok) {
            records += 1;
            bytes += rec.seq().as_ref().len() as u64;
        }
        sets.push(rs);
    }
    let produce_elapsed = produce_start.elapsed();

    if records == 0 || bytes == 0 {
        return Ok(None);
    }

    // Consumer: map in place, producer finished.
    let consume_start = Instant::now();
    for rs in &sets {
        for rec in rs.iter().filter_map(Result::ok) {
            map_one(rec.seq().as_ref());
        }
    }
    let consume_elapsed = consume_start.elapsed();

    let gb = bytes as f64 / 1e9;
    Ok(Some(InputRates {
        decode_gbps_per_stream: safe_rate(gb, produce_elapsed.as_secs_f64()),
        map_gbps_per_thread: safe_rate(gb, consume_elapsed.as_secs_f64()),
        records_sampled: records,
        bytes_sampled: bytes,
    }))
}

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
    fn preference_parses_its_command_line_forms() {
        use DecoderPreference as P;
        assert_eq!(P::parse("auto").unwrap(), P::Auto);
        assert_eq!(P::parse("serial").unwrap(), P::Serial);
        assert_eq!(P::parse("PARALLEL").unwrap(), P::Parallel { workers_per_file: None });
        assert_eq!(
            P::parse("parallel=8").unwrap(),
            P::Parallel { workers_per_file: Some(8) }
        );
        // A zero worker request is a parallel request, not a serial one.
        assert_eq!(
            P::parse("parallel=0").unwrap(),
            P::Parallel { workers_per_file: Some(1) }
        );
        assert!(P::parse("sometimes").is_err());
        assert!(P::parse("parallel=lots").is_err());
    }

    #[test]
    fn a_preference_cannot_make_a_pipe_seekable() {
        // Honouring `parallel` on a FIFO would spend coordinator threads to get
        // sequential decoding, so the forcing still wins.
        let d = preference_choice(
            DecoderPreference::Parallel { workers_per_file: Some(8) },
            InputKind::Stream,
        )
        .unwrap();
        assert!(!d.parallel);
        assert_eq!(d.reason, Reason::NonSeekableInput);

        // On a regular file it is honoured.
        let d = preference_choice(
            DecoderPreference::Parallel { workers_per_file: None },
            InputKind::Regular,
        )
        .unwrap();
        assert!(d.parallel);
        assert_eq!(d.reason, Reason::UserRequested);
    }

    #[test]
    fn auto_defers_to_the_normal_path() {
        assert!(preference_choice(DecoderPreference::Auto, InputKind::Regular).is_none());
        assert!(preference_choice(DecoderPreference::Auto, InputKind::Stream).is_none());
    }

    /// Decode and mapping must be timed separately, each with the other idle.
    /// Interleaving them measured "decode while mapping" and understated stream
    /// throughput badly.
    #[test]
    fn decode_and_map_rates_are_measured_independently() {
        let mut n = 0;
        let rates = probe(
            80,
            move || {
                n += 1;
                if n > 80 {
                    return None;
                }
                // A pathologically slow cold open, as a real one is.
                if n <= 10 {
                    std::thread::sleep(std::time::Duration::from_millis(2));
                }
                Some(vec![0u8; 1_000_000])
            },
            |_| std::thread::sleep(std::time::Duration::from_micros(50)),
        )
        .unwrap();
        // 80/8 = 10 warm-up records, so 70 are measured.
        assert_eq!(rates.records_sampled, 70);
        assert_eq!(rates.bytes_sampled, 70_000_000);
        // Decode saw no sleeps at all once warm, so it must far outpace mapping.
        assert!(
            rates.decode_gbps_per_stream > rates.map_gbps_per_thread,
            "decode {} map {}",
            rates.decode_gbps_per_stream,
            rates.map_gbps_per_thread
        );
    }

    #[test]
    fn probe_reports_none_on_empty_input() {
        assert!(probe(100, || None, |_| {}).is_none());
    }

    #[test]
    fn conflicts_are_reported_only_when_real() {
        use DecoderPreference as P;
        // Auto never conflicts, whatever the input.
        assert!(preference_conflict(P::Auto, InputKind::Stream, InputCompression::NotGzip).is_none());
        // Parallel on a pipe.
        assert_eq!(
            preference_conflict(
                P::Parallel { workers_per_file: None },
                InputKind::Stream,
                InputCompression::Unknown
            ),
            Some(PreferenceConflict::ParallelOnStream)
        );
        // Any explicit decoder on non-gzip input, including `serial`, since the
        // flag implies a choice that does not exist there.
        assert_eq!(
            preference_conflict(P::Serial, InputKind::Regular, InputCompression::NotGzip),
            Some(PreferenceConflict::NotGzipInput)
        );
        // Honourable request: nothing to say.
        assert!(
            preference_conflict(
                P::Parallel { workers_per_file: Some(4) },
                InputKind::Regular,
                InputCompression::Gzip
            )
            .is_none()
        );
    }

}
