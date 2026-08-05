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
///
/// A threads-per-file bound used to live here too, justified as a forcing on
/// the grounds that serial supply exceeded what so few threads could consume
/// "at the fastest per-thread rate ever measured". That reasoning embeds a
/// measured consumer rate, so it was empirical rather than logical, and it
/// inherited every objection to fitting a constant on one machine — a consumer
/// lighter than anything we had measured would break it. Consumer cost is a
/// continuous quantity the probe already measures, through the very closure the
/// caller supplies, so it is measured rather than assumed.
///
/// A user's explicit `--decoder serial` is also honoured unconditionally, but
/// earlier: see [`preference_choice`], which runs before this.
pub fn forced_choice(kind: InputKind, _num_files: usize, map_threads: usize) -> Option<Decision> {
    if kind == InputKind::Stream {
        return Some(Decision { parallel: false, reason: Reason::NonSeekableInput });
    }
    if map_threads <= 1 {
        return Some(Decision { parallel: false, reason: Reason::SingleThreaded });
    }
    None
}

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
        next_record()?;
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

/// Direct measurement of whether decode can keep the mapping threads fed.
///
/// # Why this replaces comparing two rates
///
/// The rate comparison asks whether `files x supply` exceeds
/// `map_threads x demand`, with both terms measured in isolation. Its supply
/// term is wrong in a way no calibration fixes: the producer never gets a whole
/// core. `paraseq` serializes `fill` behind a per-file mutex, the thread
/// holding it also maps, and every thread competes for memory bandwidth.
/// Achieved supply therefore *falls* as `-t` rises — exactly where the parallel
/// decoder earns its keep. Validated against a measured crossover surface, that
/// model got 4 of 8 points wrong, all predicting serial where parallel won.
///
/// It is also machine-specific in a way a constant cannot repair. The
/// threads-per-file constant that outperformed it was fit on one CPU, one disk
/// and one memory system; supply and demand move independently across hardware,
/// so nothing about that number transfers.
///
/// This instead runs the real shape in miniature — one producer feeding a
/// bounded queue, `map_threads` consumers draining it — and observes **which
/// side blocks**. That is the question the decision actually turns on, it is
/// measured on the hardware in front of us, and it is immune to error in either
/// absolute rate.
#[derive(Debug, Clone, Copy)]
pub struct Starvation {
    /// Fraction of consumer time spent waiting for the producer, in `0.0..=1.0`.
    ///
    /// High means the mapping threads are decode-starved and a parallel decoder
    /// has something to give. Low means decode already keeps up.
    pub consumer_wait_fraction: f64,
    /// Records mapped during the observation.
    pub records: usize,
    /// How long the observation ran.
    pub elapsed: std::time::Duration,
    /// Sampling windows observed.
    pub windows: usize,
    /// Why sampling stopped: `decisive`, `stable` or `ceiling`.
    pub stopped_because: &'static str,
}

/// Validated against the crossover surface that broke the rate model, with
/// consumers given the share of mapping threads *one file* must feed
/// (`map_threads / num_files`), since a real run fills every file concurrently:
///
/// | threads/file | consumer wait | measured parallel speedup |
/// |---:|---:|---|
/// | 4 | 1.7-1.9% | 1.09x, 1.07x, 0.92x |
/// | 8 | 44.4-46.2% | 1.92x, 1.54x, 1.20x |
/// | 16 | 57.2-58.6% | 3.05x, 2.08x |
///
/// The wait fraction is monotonic in threads-per-file and consistent *across
/// different file counts*, which is what makes it a measurement rather than a
/// restatement of the input. Six of eight points come out right; both misses
/// are in the 4 threads/file band, whose own outcomes average 1.03x. No
/// threshold can rank within that band, and none needs to.
///
/// A fitted threads-per-file constant scored 7 of 8 on the same surface, but
/// only by winning those two near-ties, and only on the machine it was fit to:
/// supply and demand move independently across CPUs, disks and memory systems.
/// The starvation fraction is measured wherever the code runs, which is why it
/// replaced the constant outright rather than supplementing it.
///
/// Consumer-wait fraction above which the parallel decoder is worth its cost.
///
/// Deliberately low. A quarter of mapping capacity lost to waiting is already
/// worth recovering, and the asymmetry established elsewhere — enabling the
/// parallel decoder needlessly costs a few percent, failing to enable it costs
/// up to 3x — argues for erring toward it.
pub const STARVATION_THRESHOLD: f64 = 0.25;

/// How long the probe samples, and when it may stop early.
///
/// The probe used to run for a fixed 150 ms. That is wrong in both directions:
/// the queue starts empty, so consumers wait at the beginning *whichever* side
/// is really bound, and a short window lets that startup transient dominate —
/// biasing the answer toward "parallel" exactly when a heavy consumer makes the
/// transient longest. Meanwhile a workload whose answer is obvious in 50 ms
/// paid the full budget anyway.
///
/// So sample in short windows and stop when the answer stops moving, or when it
/// is unambiguous, whichever comes first. Clear-cut inputs pay a fraction of the
/// ceiling; genuinely marginal ones pay more, which is where the extra
/// measurement is worth something.
#[derive(Debug, Clone, Copy)]
pub struct ProbeConfig {
    /// Length of one sampling window. Wait fractions are computed from
    /// per-window deltas, never cumulatively — a cumulative average drags the
    /// startup transient through the whole measurement.
    pub window: std::time::Duration,
    /// Windows that must elapse before any verdict, so the probe cannot
    /// converge *on* the transient it is trying to discard.
    pub min_windows: usize,
    /// Recent windows considered when testing for stability and for an
    /// unambiguous verdict.
    pub smoothing_windows: usize,
    /// Spread across the recent windows below which the answer counts as
    /// settled.
    pub tolerance: f64,
    /// Distance from [`STARVATION_THRESHOLD`] beyond which the verdict is
    /// unambiguous and further sampling cannot change it.
    pub decisive_margin: f64,
    /// Hard ceiling on total probe time.
    pub ceiling: std::time::Duration,
}

impl Default for ProbeConfig {
    fn default() -> Self {
        Self {
            window: std::time::Duration::from_millis(25),
            min_windows: 3,
            smoothing_windows: 3,
            tolerance: 0.05,
            decisive_margin: 0.20,
            ceiling: std::time::Duration::from_millis(
                env_u64("PISCEM_PROBE_CEILING_MS", 1_000),
            ),
        }
    }
}

impl ProbeConfig {
    /// Override the ceiling, in milliseconds.
    pub fn with_ceiling_ms(mut self, ms: u64) -> Self {
        self.ceiling = std::time::Duration::from_millis(ms);
        self
    }

    /// Override the sampling window, in milliseconds.
    pub fn with_window_ms(mut self, ms: u64) -> Self {
        self.window = std::time::Duration::from_millis(ms.max(1));
        self
    }
}

fn env_u64(key: &str, default: u64) -> u64 {
    std::env::var(key)
        .ok()
        .and_then(|v| v.parse().ok())
        .filter(|&v| v > 0)
        .unwrap_or(default)
}


/// Run one producer against `consumers` mapping threads and report how much of
/// the consumers' time went to waiting.
///
/// `map_one` must be callable from several threads at once; it stands in for
/// whatever the caller means by mapping one record, and *is* the consumer cost
/// — a caller with two mapping modes passes the closure for the mode it is
/// about to run, and the measurement follows continuously. There is no need to
/// classify the workload.
pub fn probe_starvation<M>(
    path: &Path,
    consumers: usize,
    cfg: ProbeConfig,
    map_one: M,
) -> anyhow::Result<Option<Starvation>>
where
    M: Fn(&[u8]) + Sync,
{
    use paraseq::fastx;
    use paraseq::prelude::*;
    use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
    use std::time::Instant;

    let (reader, _fmt) = niffler::send::from_path(path)?;
    let mut reader = fastx::Reader::new(reader)?;

    // Bounded so a fast producer blocks rather than buffering the file: the
    // queue depth is what makes "which side blocks" a meaningful question.
    let (tx, rx) = crossbeam::channel::bounded::<fastx::RecordSet>(consumers.max(1) * 2);

    let done = AtomicBool::new(false);
    let waited_nanos = AtomicU64::new(0);
    let worked_nanos = AtomicU64::new(0);
    let records = AtomicU64::new(0);

    let start = Instant::now();
    let mut windows: Vec<f64> = Vec::new();
    let mut stopped_because = "ceiling";

    std::thread::scope(|scope| {
        for _ in 0..consumers.max(1) {
            let rx = rx.clone();
            let (done, waited, worked, records) = (&done, &waited_nanos, &worked_nanos, &records);
            let map_one = &map_one;
            scope.spawn(move || {
                loop {
                    let t = Instant::now();
                    let Ok(rs) = rx.recv() else { break };
                    waited.fetch_add(t.elapsed().as_nanos() as u64, Ordering::Relaxed);

                    let t = Instant::now();
                    let mut n = 0u64;
                    for rec in rs.iter().filter_map(Result::ok) {
                        map_one(rec.seq().as_ref());
                        n += 1;
                    }
                    worked.fetch_add(t.elapsed().as_nanos() as u64, Ordering::Relaxed);
                    records.fetch_add(n, Ordering::Relaxed);
                    if done.load(Ordering::Relaxed) {
                        break;
                    }
                }
            });
        }
        drop(rx);

        // Producer on its own thread, leaving this one free to sample.
        let producer = {
            let done = &done;
            let tx = tx.clone();
            scope.spawn(move || {
                while !done.load(Ordering::Relaxed) {
                    let mut rs = reader.new_record_set();
                    match rs.fill(&mut reader) {
                        Ok(true) => {}
                        _ => break,
                    }
                    if tx.send(rs).is_err() {
                        break;
                    }
                }
            })
        };
        drop(tx);

        // Sample in windows, judging from per-window deltas. A running mean
        // would carry the startup transient -- consumers waiting on an empty
        // queue -- through the entire measurement.
        let (mut prev_wait, mut prev_work) = (0u64, 0u64);
        loop {
            std::thread::sleep(cfg.window);

            let w = waited_nanos.load(Ordering::Relaxed);
            let k = worked_nanos.load(Ordering::Relaxed);
            let (dw, dk) = (w.saturating_sub(prev_wait), k.saturating_sub(prev_work));
            prev_wait = w;
            prev_work = k;
            if dw + dk > 0 {
                windows.push(dw as f64 / (dw + dk) as f64);
            }

            if start.elapsed() >= cfg.ceiling || producer.is_finished() {
                break;
            }
            if windows.len() < cfg.min_windows.max(cfg.smoothing_windows) {
                continue;
            }

            let recent = &windows[windows.len() - cfg.smoothing_windows..];
            let mean = recent.iter().sum::<f64>() / recent.len() as f64;
            let spread = recent.iter().cloned().fold(f64::MIN, f64::max)
                - recent.iter().cloned().fold(f64::MAX, f64::min);

            // Unambiguous: no further sampling can move it across the line.
            if (mean - STARVATION_THRESHOLD).abs() > cfg.decisive_margin {
                stopped_because = "decisive";
                break;
            }
            // Settled: the answer has stopped moving.
            if spread <= cfg.tolerance {
                stopped_because = "stable";
                break;
            }
        }

        done.store(true, Ordering::Relaxed);
    });

    let n = records.load(Ordering::Relaxed) as usize;
    if n == 0 || windows.is_empty() {
        return Ok(None);
    }

    // Report the smoothed tail, so the warm-up windows -- whose only content is
    // the queue filling from empty -- do not colour the verdict.
    let tail = &windows[windows.len().saturating_sub(cfg.smoothing_windows)..];
    let fraction = tail.iter().sum::<f64>() / tail.len() as f64;

    Ok(Some(Starvation {
        consumer_wait_fraction: fraction,
        records: n,
        elapsed: start.elapsed(),
        windows: windows.len(),
        stopped_because,
    }))
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

/// How many times more records the producer side reads than the consumer maps.
///
/// The producer is fast, so a consumer-sized sample times it over ~20 ms, which
/// on a validation run read 0.83 GB/s against 0.43-0.52 measured over a whole
/// file. Ten times the records costs a few extra milliseconds and removes that.
const PRODUCER_SAMPLE_FACTOR: usize = 10;

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
    //
    // Measured over more records than the consumer needs. The two rates want
    // different sample sizes: mapping is slow enough that a few thousand
    // records take tens of milliseconds, while producing the same few thousand
    // takes ~20 ms -- short enough that the result was dominated by noise and
    // read 0.83 GB/s where a whole-file measurement of the identical work gives
    // 0.43-0.52. Reading further costs little and stabilises it.
    let produce_sample = sample.saturating_mul(PRODUCER_SAMPLE_FACTOR);
    let mut sets: Vec<fastx::RecordSet> = Vec::new();
    let mut retained_records = 0usize;
    let mut records = 0usize;
    let mut bytes = 0u64;
    let produce_start = Instant::now();
    while records < produce_sample {
        let mut rs = reader.new_record_set();
        match rs.fill(&mut reader) {
            Ok(true) => {}
            _ => break,
        }
        for rec in rs.iter().filter_map(Result::ok) {
            records += 1;
            bytes += rec.seq().as_ref().len() as u64;
        }
        // Retain only what the consumer phase will map; drop the rest as we go.
        // Holding all of them grows the working set through the timed loop and
        // cost ~20% of the measured producer rate -- 0.70 GB/s where dropping
        // gives 0.88, against a whole-file truth of 0.82.
        if retained_records < sample {
            retained_records += rs.iter().filter_map(Result::ok).count();
            sets.push(rs);
        }
    }
    let produce_elapsed = produce_start.elapsed();

    if records == 0 || bytes == 0 {
        return Ok(None);
    }

    // Consumer: map in place, producer finished. Only the first `sample`
    // records -- mapping is the expensive half, and its rate is stable well
    // before the producer's is.
    let mut mapped = 0usize;
    let mut mapped_bytes = 0u64;
    let consume_start = Instant::now();
    'consume: for rs in &sets {
        for rec in rs.iter().filter_map(Result::ok) {
            if mapped >= sample {
                break 'consume;
            }
            let seq = rec.seq();
            mapped_bytes += seq.as_ref().len() as u64;
            map_one(seq.as_ref());
            mapped += 1;
        }
    }
    let consume_elapsed = consume_start.elapsed();

    Ok(Some(InputRates {
        decode_gbps_per_stream: safe_rate(bytes as f64 / 1e9, produce_elapsed.as_secs_f64()),
        map_gbps_per_thread: safe_rate(
            mapped_bytes as f64 / 1e9,
            consume_elapsed.as_secs_f64(),
        ),
        records_sampled: mapped,
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

    /// Only genuinely logical cases are forced. A threads-per-file bound used
    /// to live here; it encoded consumer cost as a constant, which the probe
    /// measures directly, so 8 files at 16 threads must now be *measured*
    /// rather than assumed.
    #[test]
    fn thread_ratios_are_measured_not_forced() {
        assert!(forced_choice(InputKind::Regular, 8, 16).is_none());
        assert!(forced_choice(InputKind::Regular, 2, 4).is_none());
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
