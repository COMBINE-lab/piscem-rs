//! FASTX reader helpers — open and concatenate FASTQ files with decompression.
//!
//! Provides `open_concatenated_readers()` for creating a single `Read` stream
//! from multiple (possibly compressed) FASTQ files, suitable for feeding to
//! paraseq's `Reader`.

use anyhow::{Context, Result};
use std::path::Path;

// Re-export paraseq types used by the mapping pipeline.
pub use paraseq::Record;
pub use paraseq::fastq;
pub use paraseq::fastx::{Collection, CollectionType};
pub use paraseq::parallel::{
    MultiParallelProcessor, PairedParallelProcessor, ParallelProcessor, ParallelReader,
};

/// Batch size for paraseq record sets. paraseq's default is 1024; we raise it
/// to reduce contention on the per-reader mutex when many mapping threads drain
/// batches faster than a single thread can fill them.
///
/// Swept and confirmed on the Flex pair at `-t 32` using the real pipeline at
/// three consumer costs. 4096 and 16384 form a
/// plateau, differing by -0.8% to +3.8% — inside the ~3% noise floor, so neither
/// is preferable. Outside it: **1024 is 8-10% worse** and **65536 is 4-9%
/// worse**. Starvation also *falls* as batches grow (55.5% -> 41.6% at the
/// marginal consumer cost), so the mutex is better amortised by larger batches
/// rather than blocked behind them.
///
/// A short measurement says the opposite, and convincingly. At a 75 ms probe,
/// throughput fell monotonically with batch size and starvation rose — because
/// 65536 records x 32 threads is more groups than the probe reads, so most
/// workers never started and the comparison measured startup. Sweep this only
/// at steady state.
pub const READER_BATCH_SIZE: usize = 16384;

/// scATAC batch size, chosen to preserve the resizable-pool acknowledgement
/// bound under its much more expensive per-record mapping callback.
///
/// A worker can retire only after completing the batch it owns; 16K scATAC
/// records exceeded the broker's two-second drain timeout on a real
/// negative-scaling workload. Other modalities retain [`READER_BATCH_SIZE`]
/// because their measured throughput plateau favors the larger batch.
pub const SCATAC_READER_BATCH_SIZE: usize = 1024;

/// Build a paraseq fastx Reader with the shared piscem batch size.
pub fn reader_with_batch_size<R: std::io::Read>(
    inner: R,
) -> Result<paraseq::fastx::Reader<R>, paraseq::Error> {
    paraseq::fastx::Reader::new_with_batch_size(inner, READER_BATCH_SIZE)
}

/// Build a scATAC reader with resize-bounded work granularity.
pub fn scatac_reader<R: std::io::Read>(
    inner: R,
) -> Result<paraseq::fastx::Reader<R>, paraseq::Error> {
    paraseq::fastx::Reader::new_with_batch_size(inner, SCATAC_READER_BATCH_SIZE)
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Try to open `path` with rapidgzip's parallel gzip decoder.
///
/// Returns `Ok(None)` when the file is not gzip, so the caller falls back to
/// niffler (which also covers zstd/bz2/xz and uncompressed input).
///
/// Motivation: paraseq holds one mutex per input file, so exactly one thread
/// can be inflating a given file at any moment — a hard ceiling of two inflate
/// streams for a paired run, regardless of `-t`. On the full Flex dataset that
/// caps the run at ~11 min while the mapping threads sit ~80% idle. rapidgzip
/// decodes a *single* gzip member in parallel internally, so the ceiling moves
/// without changing paraseq's structure at all.
/// An opened input plus, when the parallel gzip decoder is in use, the control
/// handle for it.
///
/// The handle must be taken *before* the reader is moved into paraseq (which
/// consumes it as `impl Read`), so it is surfaced here rather than fetched later.
pub struct OpenedInput {
    pub reader: Box<dyn std::io::Read + Send>,
    #[cfg(feature = "rapidgzip")]
    pub handle: Option<rapidgzip_core::DecoderHandle>,
}

#[cfg(feature = "rapidgzip")]
fn open_gz_rapidgzip(path: &Path, ceiling: usize) -> Result<Option<OpenedInput>> {
    use std::io::Read;

    // Sniff the magic rather than trusting the extension.
    let mut probe =
        std::fs::File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    let mut magic = [0u8; 2];
    if probe.read_exact(&mut magic).is_err() || magic != [0x1f, 0x8b] {
        return Ok(None);
    }
    drop(probe);

    // `decoder_threads` is an immutable *ceiling*, not an allocation: workers
    // are created lazily as the decoder's own calibrator finds useful work, and
    // retire again under sustained consumer backpressure. So this can be set
    // generously; `DecodeBudget` caps what the process actually spends.
    let decoder = rapidgzip_core::Decoder::builder()
        .decoder_threads(ceiling)
        .build()
        .map_err(|e| anyhow::anyhow!("rapidgzip-rust decoder config: {e}"))?;
    let reader = decoder
        .open(path)
        .map_err(|e| anyhow::anyhow!("rapidgzip-rust failed to open {}: {e}", path.display()))?;
    let handle = reader.handle();

    // Deliberately do NOT call `set_worker_limit` here.
    //
    // On single-member gzip the application limit is not a ceiling the decoder
    // can grow past later — it is a one-shot path selector, read once while
    // choosing between the parallel and sequential backends. `should_probe`
    // requires at least two effective workers; below that rapidgzip commits to
    // `DecoderPath::Sequential`, pins its adaptive target to 1, and never
    // spawns or registers a worker for the rest of the file. Raising the limit
    // afterwards is a no-op.
    //
    // A previous version set an "initial limit" of 1 here, believing it merely
    // a cheap starting point a controller could ratchet up. It was instead a
    // latch that disabled parallel decoding outright: the decoder reported
    // `Idle`/`spawned_workers = 0` for entire runs while the mapping threads
    // sat 77% starved. `RuntimeState::new` already seeds the runtime limit from
    // `decoder_threads`, so leaving it alone is both correct and simpler.
    //
    // There is also a race that makes any open-time clamp unsound: `open`
    // spawns the coordinator, which may begin classifying before this line
    // could run.
    tracing::debug!(
        "rapidgzip-rust: {} (ceiling {} workers)",
        path.display(),
        ceiling
    );
    Ok(Some(OpenedInput {
        reader: Box::new(reader),
        handle: Some(handle),
    }))
}

/// Open a single file with automatic decompression, surfacing the decoder
/// control handle when the parallel gzip path is taken.
pub fn open_input(path: impl AsRef<Path>, ceiling: usize) -> Result<OpenedInput> {
    let path = path.as_ref();

    // `ceiling == 0` selects the serial path explicitly, independent of whether
    // the feature is on.
    //
    // Non-regular inputs are excluded here rather than left to the decoder.
    // `rapidgzip` gates its parallel path on `file_type().is_file()` and falls
    // back to sequential decoding otherwise, so it has nothing to offer on a
    // FIFO -- and reaching it is actively harmful, because the magic sniff
    // below opens the path, consumes two bytes, and closes it, after which
    // re-opening a FIFO blocks forever waiting for a writer that has gone.
    // Process substitution (`-r <(zcat ...)`) hits exactly this.
    #[cfg(feature = "rapidgzip")]
    if ceiling > 0
        && crate::io::calibrate::classify_input(path) == crate::io::calibrate::InputKind::Regular
        && let Some(opened) = open_gz_rapidgzip(path, ceiling)?
    {
        return Ok(opened);
    }
    #[cfg(not(feature = "rapidgzip"))]
    let _ = ceiling;

    let (reader, _format) = niffler::send::from_path(path)
        .with_context(|| format!("failed to open {}", path.display()))?;
    Ok(OpenedInput {
        reader,
        #[cfg(feature = "rapidgzip")]
        handle: None,
    })
}

/// Open a gzip file attached to a shared decode pool.
///
/// Every input in a run shares one pool, so the decode budget is a single
/// aggregate rather than a per-file share that has to be divided up front. That
/// division is what produced every thread-accounting bug in this area: a
/// decoder granted more than it uses starves its peers, and one granted less
/// than it needs cannot borrow.
///
/// `decoder_threads` is per-operation *headroom*, not an allocation. It bounds
/// what one decoder may request and is what admission reads when choosing a
/// backend, so it should be generous — the pool's mutable limit is what
/// actually constrains concurrency, and it no longer affects path selection.
///
/// Returns `Ok(None)` when the file is not gzip, so the caller falls back to
/// niffler.
#[cfg(feature = "rapidgzip")]
pub fn open_gz_pooled(
    path: &Path,
    pool: &rapidgzip_core::DecoderPool,
    decoder_threads: usize,
) -> Result<Option<OpenedInput>> {
    use std::io::Read;

    let mut probe =
        std::fs::File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    let mut magic = [0u8; 2];
    if probe.read_exact(&mut magic).is_err() || magic != [0x1f, 0x8b] {
        return Ok(None);
    }
    drop(probe);

    let decoder = rapidgzip_core::Decoder::builder()
        .decoder_threads(decoder_threads.max(2))
        .decoder_pool(pool.clone())
        .build()
        .map_err(|e| anyhow::anyhow!("rapidgzip-rust decoder config: {e}"))?;
    let reader = decoder
        .open(path)
        .map_err(|e| anyhow::anyhow!("rapidgzip-rust failed to open {}: {e}", path.display()))?;
    let handle = reader.handle();
    tracing::debug!(
        "rapidgzip-rust (pooled): {} (headroom {} workers)",
        path.display(),
        decoder_threads
    );
    Ok(Some(OpenedInput {
        reader: Box::new(reader),
        handle: Some(handle),
    }))
}

/// Open an input for a pooled run, falling back to the serial path where the
/// parallel decoder cannot help.
///
/// The per-file downgrade lives here: a FIFO among the inputs does not cost the
/// regular files their parallel decoder, and a non-gzip input never had a
/// choice to make.
#[cfg(feature = "rapidgzip")]
pub fn open_input_pooled(
    path: impl AsRef<Path>,
    pool: Option<&rapidgzip_core::DecoderPool>,
    decoder_threads: usize,
) -> Result<OpenedInput> {
    let path = path.as_ref();
    if let Some(pool) = pool
        && crate::io::calibrate::classify_input(path) == crate::io::calibrate::InputKind::Regular
        && let Some(opened) = open_gz_pooled(path, pool, decoder_threads)?
    {
        return Ok(opened);
    }
    open_input(path, 0)
}

/// Open a single file with automatic decompression (gzip, zstd, etc.).
///
/// Convenience wrapper for call sites that do not participate in decoder
/// budgeting (bulk and scATAC); it discards the control handle.
pub(crate) fn open_with_decompression(
    path: impl AsRef<Path>,
) -> Result<Box<dyn std::io::Read + Send>> {
    let c = default_decoder_ceiling();
    Ok(open_input(path, c)?.reader)
}

/// Per-file decoder ceiling when the caller has no budget of its own.
pub(crate) fn default_decoder_ceiling() -> usize {
    std::env::var("PISCEM_RAPIDGZIP_THREADS")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(16)
}

/// Where a fixed decoder allocation came from.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize)]
#[serde(rename_all = "snake_case")]
pub enum DecodePinSource {
    /// `--decoder parallel=N`, where `N` is a per-file request.
    CliPerFile,
    /// `PISCEM_DECODE_SLOTS=N`, an aggregate measurement/reproduction control.
    AggregateEnvironment,
    /// The legacy `PISCEM_RAPIDGZIP_THREADS=N` per-file control.
    LegacyPerFileEnvironment,
}

/// How the run will allocate decompression slots.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum DecodeAllocation {
    /// No parallel decoder pool; mapping receives the entire effective budget.
    Serial,
    /// The thread broker may move slots between decoding and mapping.
    Adaptive,
    /// A fixed aggregate number of decode slots.
    PinnedAggregate {
        slots: usize,
        source: DecodePinSource,
    },
    /// A fixed request per decoder-capable input. Reconciled against the number
    /// of decoder handles actually opened before mapping starts.
    PinnedPerFile {
        slots_per_file: usize,
        source: DecodePinSource,
    },
}

/// One authoritative plan for every thread-related component of a mapping run.
///
/// Pool maxima, opening targets, the broker budget, logs, and telemetry must all
/// use `effective_budget`; retaining the user's original request separately is
/// what makes a cgroup/cpuset reduction visible without accidentally undoing it.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize)]
pub struct ExecutionPlan {
    /// The command-line budget before resource-environment capping.
    pub requested_budget: usize,
    /// The budget actually divided by this run.
    pub effective_budget: usize,
    /// Threads that will map. With `decode_slots`, this equals
    /// `effective_budget` in every non-serial plan.
    pub map_threads: usize,
    /// Aggregate decoder execution slots.
    pub decode_slots: usize,
    /// Requested decoder allocation and its pin semantics.
    pub allocation: DecodeAllocation,
    /// Decode slots requested before a fixed request was clamped to leave one
    /// mapping thread. `None` for serial and adaptive allocation.
    pub requested_decode_slots: Option<usize>,
}

/// Thread-planning and adaptive-control telemetry returned by a mapping
/// pipeline for inclusion in `map_info.json`.
#[derive(Debug)]
pub struct PipelineOutcome {
    pub execution_plan: ExecutionPlan,
    pub broker_report: Option<thread_broker::BrokerReport>,
    pub producer_measurement: Option<thread_broker::ProducerMeasurementStats>,
    pub consumer_measurement: crate::io::threads::ConsumerMeasurement,
    /// Real mapping-pipeline wall time, excluding index load and output
    /// backpatching, for reproducible fixed-split comparisons.
    pub mapping_elapsed_secs: f64,
}

/// Controller cadence chosen from the effective execution budget.
///
/// Small-budget jobs are often short enough to end during the general 400 ms
/// warm-up. The observed consumer startup transient clears at about 75 ms, so a
/// 100 ms warm-up followed by three 25 ms clean windows lets `-t <= 8` adapt
/// before a quarter-second job finishes. Larger runs keep the conservative
/// defaults used for the original convergence measurements.
pub fn broker_config_for_budget(effective_budget: usize) -> thread_broker::BrokerConfig {
    if effective_budget <= 8 {
        thread_broker::BrokerConfig {
            sample_interval: std::time::Duration::from_millis(25),
            warmup: std::time::Duration::from_millis(100),
            ..thread_broker::BrokerConfig::default()
        }
    } else {
        thread_broker::BrokerConfig::default()
    }
}

impl ExecutionPlan {
    /// Whether a shared parallel decoder pool should be constructed.
    pub fn parallel_gzip(self) -> bool {
        !matches!(self.allocation, DecodeAllocation::Serial)
    }

    /// Whether the broker, rather than the user, owns the split.
    pub fn adaptive(self) -> bool {
        matches!(self.allocation, DecodeAllocation::Adaptive)
    }

    /// Recompute the final allocation from the decoder handles actually opened.
    ///
    /// This returns reserved mapping threads when regular inputs turn out not to
    /// be gzip, and makes `parallel=N` genuinely per decoder-capable file rather
    /// than per path listed on the command line.
    pub fn reconcile_parallel_decoders(&mut self, decoders: usize) {
        if decoders == 0 {
            self.allocation = DecodeAllocation::Serial;
            self.decode_slots = 0;
            self.map_threads = self.effective_budget;
            self.requested_decode_slots = None;
            return;
        }

        if let DecodeAllocation::PinnedPerFile { slots_per_file, .. } = self.allocation {
            let requested = slots_per_file.saturating_mul(decoders);
            self.set_fixed_decode_slots(requested);
        }
    }

    fn set_fixed_decode_slots(&mut self, requested: usize) {
        debug_assert!(self.effective_budget >= 2);
        let maximum = self.effective_budget - 1;
        let applied = requested.clamp(1, maximum);
        self.requested_decode_slots = Some(requested);
        self.decode_slots = applied;
        self.map_threads = self.effective_budget - applied;
    }
}

/// A coarse opening split of `total_threads` between mapping and decoding.
///
/// # This is a starting point, not an answer
///
/// It used to be the answer, computed from rates measured before the run. That
/// is gone: measured against the optimum on our own workloads, the closed-form
/// model was **60% off** and a two-point measurement 7-38% off depending on
/// where the second point was taken, with no placement good on both. A plain
/// constant beat both.
///
/// Both sides are now resizable mid-run, so the split is converged to instead —
/// see `thread_broker`. All this has to do is start somewhere sane and not
/// waste the first few sampling intervals.
///
/// **Biased toward mapping on purpose.** Growth is cheap and safe now, so
/// starting low costs only convergence time, while starting high wastes threads
/// outright on the many workloads where decode never binds — measured at 10.0%
/// consumer idle on paired bulk RNA-seq, where the right decode allocation is
/// nearly none.
pub fn plan_thread_budget(
    total_threads: usize,
    num_files: usize,
    parallel_selected: bool,
    decoder_preference: crate::io::calibrate::DecoderPreference,
) -> Result<ExecutionPlan> {
    let requested = total_threads.max(1);
    // `-t` is a request, not a guarantee: under a cpuset or a cgroup CPU quota
    // the process may hold far fewer cores than it names, and budgeting against
    // the larger number just oversubscribes. (Bound suggested by salmon #1072.)
    let usable = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(requested);

    let allocation = requested_allocation(parallel_selected, decoder_preference)?;
    let plan = build_execution_plan(requested, usable, num_files, allocation);
    if let Some(requested_slots) = plan.requested_decode_slots
        && requested_slots != plan.decode_slots
    {
        tracing::warn!(
            "requested {requested_slots} aggregate decoder slots, but the effective \
             budget {} can provide at most {}; using {} decode / {} mapping",
            plan.effective_budget,
            plan.effective_budget.saturating_sub(1),
            plan.decode_slots,
            plan.map_threads,
        );
    }
    Ok(plan)
}

fn requested_allocation(
    parallel_selected: bool,
    decoder_preference: crate::io::calibrate::DecoderPreference,
) -> Result<DecodeAllocation> {
    if !cfg!(feature = "rapidgzip") || !parallel_selected {
        return Ok(DecodeAllocation::Serial);
    }

    if let crate::io::calibrate::DecoderPreference::Parallel {
        workers_per_file: Some(slots_per_file),
    } = decoder_preference
    {
        if slots_per_file == 0 {
            anyhow::bail!("--decoder parallel=N requires N to be at least 1");
        }
        return Ok(DecodeAllocation::PinnedPerFile {
            slots_per_file,
            source: DecodePinSource::CliPerFile,
        });
    }

    let aggregate = std::env::var("PISCEM_DECODE_SLOTS").ok();
    let legacy = std::env::var("PISCEM_RAPIDGZIP_THREADS").ok();
    if let Some(allocation) = environment_allocation(aggregate.as_deref(), legacy.as_deref())? {
        return Ok(allocation);
    }

    Ok(DecodeAllocation::Adaptive)
}

fn environment_allocation(
    aggregate: Option<&str>,
    legacy: Option<&str>,
) -> Result<Option<DecodeAllocation>> {
    if aggregate.is_some() && legacy.is_some() {
        anyhow::bail!(
            "PISCEM_DECODE_SLOTS and PISCEM_RAPIDGZIP_THREADS are both set; \
             choose the aggregate or per-file control, not both"
        );
    }
    if let Some(raw) = aggregate {
        let slots = raw.parse::<usize>().map_err(|_| {
            anyhow::anyhow!("PISCEM_DECODE_SLOTS must be a positive integer, got `{raw}`")
        })?;
        if slots == 0 {
            anyhow::bail!("PISCEM_DECODE_SLOTS must be at least 1");
        }
        return Ok(Some(DecodeAllocation::PinnedAggregate {
            slots,
            source: DecodePinSource::AggregateEnvironment,
        }));
    }
    if let Some(raw) = legacy {
        let slots_per_file = raw.parse::<usize>().map_err(|_| {
            anyhow::anyhow!("PISCEM_RAPIDGZIP_THREADS must be a positive integer, got `{raw}`")
        })?;
        if slots_per_file == 0 {
            anyhow::bail!("PISCEM_RAPIDGZIP_THREADS must be at least 1");
        }
        return Ok(Some(DecodeAllocation::PinnedPerFile {
            slots_per_file,
            source: DecodePinSource::LegacyPerFileEnvironment,
        }));
    }
    Ok(None)
}

fn build_execution_plan(
    requested_budget: usize,
    usable_budget: usize,
    num_files: usize,
    allocation: DecodeAllocation,
) -> ExecutionPlan {
    let requested_budget = requested_budget.max(1);
    let effective_budget = requested_budget.min(usable_budget.max(1));

    if effective_budget < 2 || matches!(allocation, DecodeAllocation::Serial) {
        return ExecutionPlan {
            requested_budget,
            effective_budget,
            map_threads: effective_budget,
            decode_slots: 0,
            allocation: DecodeAllocation::Serial,
            requested_decode_slots: None,
        };
    }

    // A two-slot opening at small budgets is the measured minimax choice: a
    // mapping-heavy short run can reclaim one slot after its first clean model,
    // while an allocation-dependent decoder may be unable to reveal that a
    // second slot is useful when observed only at the one-slot floor.
    let opening_divisor = if effective_budget <= 8 { 4 } else { 8 };
    let initial = (effective_budget / opening_divisor).clamp(1, effective_budget - 1);
    let mut plan = ExecutionPlan {
        requested_budget,
        effective_budget,
        map_threads: effective_budget - initial,
        decode_slots: initial,
        allocation,
        requested_decode_slots: None,
    };
    match allocation {
        DecodeAllocation::Serial => unreachable!("handled above"),
        DecodeAllocation::Adaptive => {}
        DecodeAllocation::PinnedAggregate { slots, .. } => plan.set_fixed_decode_slots(slots),
        DecodeAllocation::PinnedPerFile { slots_per_file, .. } => {
            let requested = slots_per_file.saturating_mul(num_files.max(1));
            plan.set_fixed_decode_slots(requested);
        }
    }
    plan
}

/// Minimum aggregate mapping target imposed by paraseq's concurrent shares.
///
/// A `Collection` fixes its batch geometry from the opening mapping target and
/// gives every concurrently active reader group a pool share. Each share has a
/// hard floor of one worker, so the broker must retain at least the maximum
/// number of shares in a batch. `arity` is one for single-end input, two for a
/// paired collection, and the modality-specific group width for multi input.
pub fn collection_share_floor(
    total_inputs: usize,
    arity: usize,
    opening_mapping_threads: usize,
) -> usize {
    let arity = arity.max(1);
    let groups = (total_inputs / arity).max(1);
    let threads = opening_mapping_threads.max(1);
    let threads_per_group = if arity == 1 {
        (threads / groups).max(1)
    } else if threads >= arity {
        (threads / groups).max(arity)
    } else {
        (threads / groups).max(1)
    };
    (threads / threads_per_group).max(1).min(groups)
}

/// Open multiple files and concatenate them into a single reader.
///
/// Automatically detects and decompresses gzip, zstd, and other formats
/// via niffler.
pub fn open_concatenated_readers(paths: &[String]) -> Result<Box<dyn std::io::Read + Send>> {
    use std::io::Read;

    if paths.is_empty() {
        anyhow::bail!("no input files specified");
    }
    if paths.len() == 1 {
        return open_with_decompression(&paths[0]);
    }
    // Collect all file readers, then chain them via a wrapper.
    let mut readers: Vec<Box<dyn Read + Send>> = Vec::with_capacity(paths.len());
    for path in paths {
        readers.push(open_with_decompression(path)?);
    }
    Ok(Box::new(MultiReader {
        readers,
        current: 0,
    }))
}

/// Concatenating reader over multiple boxed readers.
struct MultiReader {
    readers: Vec<Box<dyn std::io::Read + Send>>,
    current: usize,
}

impl std::io::Read for MultiReader {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        while self.current < self.readers.len() {
            let n = self.readers[self.current].read(buf)?;
            if n > 0 {
                return Ok(n);
            }
            self.current += 1;
        }
        Ok(0)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::{
        DecodeAllocation, DecodePinSource, build_execution_plan, collection_share_floor,
        environment_allocation,
    };

    #[test]
    fn test_open_concatenated_readers_empty() {
        let result = super::open_concatenated_readers(&[]);
        assert!(result.is_err());
    }

    #[test]
    fn serial_uses_the_whole_effective_budget() {
        let p = build_execution_plan(32, 12, 8, DecodeAllocation::Serial);
        assert_eq!(p.requested_budget, 32);
        assert_eq!(p.effective_budget, 12);
        assert_eq!(p.map_threads, 12);
        assert_eq!(p.decode_slots, 0);
        assert_eq!(p.map_threads + p.decode_slots, p.effective_budget);
    }

    #[test]
    fn adaptive_divides_the_capped_budget() {
        let p = build_execution_plan(64, 16, 2, DecodeAllocation::Adaptive);
        assert_eq!(p.effective_budget, 16);
        assert_eq!(p.map_threads, 14);
        assert_eq!(p.decode_slots, 2);
        assert_eq!(p.map_threads + p.decode_slots, 16);
        assert!(p.adaptive());
    }

    #[test]
    fn one_thread_is_always_serial() {
        for allocation in [
            DecodeAllocation::Adaptive,
            DecodeAllocation::PinnedAggregate {
                slots: 8,
                source: DecodePinSource::AggregateEnvironment,
            },
            DecodeAllocation::PinnedPerFile {
                slots_per_file: 4,
                source: DecodePinSource::CliPerFile,
            },
        ] {
            let p = build_execution_plan(1, 128, 24, allocation);
            assert_eq!(p.allocation, DecodeAllocation::Serial);
            assert_eq!((p.map_threads, p.decode_slots), (1, 0));
        }
    }

    #[test]
    fn aggregate_pin_is_clamped_but_remains_budgeted() {
        let p = build_execution_plan(
            8,
            64,
            2,
            DecodeAllocation::PinnedAggregate {
                slots: 99,
                source: DecodePinSource::AggregateEnvironment,
            },
        );
        assert_eq!(p.requested_decode_slots, Some(99));
        assert_eq!((p.map_threads, p.decode_slots), (1, 7));
        assert_eq!(p.map_threads + p.decode_slots, p.effective_budget);
    }

    #[test]
    fn per_file_pin_reconciles_against_real_decoders() {
        let mut p = build_execution_plan(
            32,
            32,
            8,
            DecodeAllocation::PinnedPerFile {
                slots_per_file: 3,
                source: DecodePinSource::CliPerFile,
            },
        );
        assert_eq!(p.requested_decode_slots, Some(24));
        p.reconcile_parallel_decoders(2);
        assert_eq!(p.requested_decode_slots, Some(6));
        assert_eq!((p.map_threads, p.decode_slots), (26, 6));
    }

    #[test]
    fn no_real_decoders_returns_every_slot_to_mapping() {
        let mut p = build_execution_plan(32, 32, 2, DecodeAllocation::Adaptive);
        p.reconcile_parallel_decoders(0);
        assert_eq!(p.allocation, DecodeAllocation::Serial);
        assert_eq!((p.map_threads, p.decode_slots), (32, 0));
    }

    #[test]
    fn environment_controls_are_unambiguous_and_validated() {
        assert!(environment_allocation(Some("3"), Some("4")).is_err());
        assert!(environment_allocation(Some("0"), None).is_err());
        assert!(environment_allocation(Some("nope"), None).is_err());
        assert_eq!(
            environment_allocation(Some("7"), None).unwrap(),
            Some(DecodeAllocation::PinnedAggregate {
                slots: 7,
                source: DecodePinSource::AggregateEnvironment,
            })
        );
        assert_eq!(
            environment_allocation(None, Some("5")).unwrap(),
            Some(DecodeAllocation::PinnedPerFile {
                slots_per_file: 5,
                source: DecodePinSource::LegacyPerFileEnvironment,
            })
        );
    }

    #[test]
    fn every_plan_respects_the_effective_budget() {
        let allocations = [
            DecodeAllocation::Serial,
            DecodeAllocation::Adaptive,
            DecodeAllocation::PinnedAggregate {
                slots: 1,
                source: DecodePinSource::AggregateEnvironment,
            },
            DecodeAllocation::PinnedAggregate {
                slots: usize::MAX,
                source: DecodePinSource::AggregateEnvironment,
            },
            DecodeAllocation::PinnedPerFile {
                slots_per_file: 3,
                source: DecodePinSource::CliPerFile,
            },
        ];

        for requested in [1, 2, 8, 32] {
            for usable in [1, 2, 6, 64] {
                for files in [0, 1, 2, 17] {
                    for allocation in allocations {
                        let plan = build_execution_plan(requested, usable, files, allocation);
                        assert_eq!(plan.requested_budget, requested);
                        assert_eq!(plan.effective_budget, requested.min(usable));
                        assert_eq!(
                            plan.map_threads + plan.decode_slots,
                            plan.effective_budget,
                            "requested={requested}, usable={usable}, files={files}, allocation={allocation:?}",
                        );
                        assert!(plan.map_threads >= 1);
                        assert!(plan.decode_slots < plan.effective_budget);
                        assert_eq!(plan.parallel_gzip(), plan.decode_slots > 0);
                    }
                }
            }
        }
    }

    #[test]
    fn execution_plan_has_stable_machine_readable_fields() {
        let plan = build_execution_plan(
            32,
            12,
            3,
            DecodeAllocation::PinnedAggregate {
                slots: 4,
                source: DecodePinSource::AggregateEnvironment,
            },
        );
        let value = serde_json::to_value(plan).unwrap();
        assert_eq!(value["requested_budget"], 32);
        assert_eq!(value["effective_budget"], 12);
        assert_eq!(value["map_threads"], 8);
        assert_eq!(value["decode_slots"], 4);
        assert_eq!(value["allocation"]["kind"], "pinned_aggregate");
        assert_eq!(value["allocation"]["slots"], 4);
        assert_eq!(value["allocation"]["source"], "aggregate_environment");
        assert_eq!(value["requested_decode_slots"], 4);
    }

    #[test]
    fn collection_floor_matches_the_number_of_concurrent_shares() {
        assert_eq!(collection_share_floor(8, 1, 28), 8);
        assert_eq!(collection_share_floor(16, 2, 12), 6);
        assert_eq!(collection_share_floor(24, 3, 10), 3);
        assert_eq!(collection_share_floor(2, 2, 1), 1);
        assert_eq!(collection_share_floor(0, 1, 8), 1);
    }

    #[test]
    fn small_budget_cadence_can_act_before_a_quarter_second() {
        let small = super::broker_config_for_budget(8);
        assert_eq!(small.sample_interval, std::time::Duration::from_millis(25));
        assert_eq!(small.warmup, std::time::Duration::from_millis(100));
        assert!(
            small.warmup + small.sample_interval * (small.smoothing_windows as u32)
                < std::time::Duration::from_millis(250)
        );

        let normal = super::broker_config_for_budget(32);
        assert_eq!(
            normal.sample_interval,
            thread_broker::BrokerConfig::default().sample_interval
        );
        assert_eq!(normal.warmup, thread_broker::BrokerConfig::default().warmup);
    }
}
