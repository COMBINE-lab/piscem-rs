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
pub const READER_BATCH_SIZE: usize = 16384;

/// Build a paraseq fastx Reader with the shared piscem batch size.
pub fn reader_with_batch_size<R: std::io::Read>(
    inner: R,
) -> Result<paraseq::fastx::Reader<R>, paraseq::Error> {
    paraseq::fastx::Reader::new_with_batch_size(inner, READER_BATCH_SIZE)
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
fn open_gz_rapidgzip(
    path: &Path,
    ceiling: usize,
    initial_limit: usize,
) -> Result<Option<OpenedInput>> {
    use std::io::Read;

    // Sniff the magic rather than trusting the extension.
    let mut probe = std::fs::File::open(path)
        .with_context(|| format!("failed to open {}", path.display()))?;
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

    // Apply the fair-share limit *here*, not once the supervisor starts.
    // `reader_with_batch_size` makes paraseq read a first batch immediately, so
    // decoding — and worker growth to the ceiling — begins before any later
    // clamp could run. Threads never created cost nothing; threads retired
    // moments after spawning have already cost their setup.
    if initial_limit > 0
        && let Err(e) = handle.set_worker_limit(initial_limit.min(ceiling))
    {
        tracing::warn!("could not set initial decoder limit: {e}");
    }
    tracing::debug!(
        "rapidgzip-rust: {} (ceiling {} workers, initial limit {})",
        path.display(),
        ceiling,
        initial_limit
    );
    Ok(Some(OpenedInput {
        reader: Box::new(reader),
        handle: Some(handle),
    }))
}

/// Open a single file with automatic decompression, surfacing the decoder
/// control handle when the parallel gzip path is taken.
pub fn open_input(
    path: impl AsRef<Path>,
    ceiling: usize,
    initial_limit: usize,
) -> Result<OpenedInput> {
    let path = path.as_ref();

    // `ceiling == 0` selects the serial path explicitly (see
    // `MIN_THREADS_PER_FILE`), independent of whether the feature is on.
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
        && let Some(opened) = open_gz_rapidgzip(path, ceiling, initial_limit)?
    {
        return Ok(opened);
    }
    #[cfg(not(feature = "rapidgzip"))]
    let _ = (ceiling, initial_limit);

    let (reader, _format) = niffler::send::from_path(path)
        .with_context(|| format!("failed to open {}", path.display()))?;
    Ok(OpenedInput {
        reader,
        #[cfg(feature = "rapidgzip")]
        handle: None,
    })
}

/// Open a single file with automatic decompression (gzip, zstd, etc.).
///
/// Convenience wrapper for call sites that do not participate in decoder
/// budgeting (bulk and scATAC); it discards the control handle.
pub(crate) fn open_with_decompression(
    path: impl AsRef<Path>,
) -> Result<Box<dyn std::io::Read + Send>> {
    let c = default_decoder_ceiling();
    Ok(open_input(path, c, c)?.reader)
}

/// Per-file decoder ceiling when the caller has no budget of its own.
pub(crate) fn default_decoder_ceiling() -> usize {
    std::env::var("PISCEM_RAPIDGZIP_THREADS")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(16)
}

fn env_usize(key: &str, default: usize) -> usize {
    std::env::var(key)
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(default)
}

/// How a run intends to spend threads on decompression.
#[derive(Debug, Clone, Copy)]
pub struct ThreadBudget {
    /// Total decoder workers all files may hold at once.
    pub decode_budget: usize,
    /// Immutable per-decoder ceiling. Generous on purpose: it is a lazy cap,
    /// and `DecodeBudget` enforces the real limit at runtime.
    pub per_file_ceiling: usize,
    /// Runtime limit each decoder starts at, so the budget holds from the very
    /// first batch rather than after the supervisor's first sample.
    pub initial_per_file: usize,
    /// Whether to use the parallel gzip decoder at all. When false, gzip input
    /// takes the serial niffler path and no decoder threads are spawned.
    pub parallel_gzip: bool,
}

/// Re-export of the single threads-per-file threshold, which lives in
/// [`crate::io::calibrate`] alongside the measurements that set it.
///
/// `plan_thread_budget` applies it as a fast path so a run that is obviously
/// serial never builds mapping state for a probe; `calibrate::forced_choice`
/// applies the same bound as its `AmpleSerialSupply` arm.
use crate::io::calibrate::MIN_THREADS_PER_FILE;

/// Decide the decompression share of the thread budget.
///
/// Measured on a 2.30 B read-pair archive (single .gz pair, ~186 KB gzip
/// members, 2 x 149.5 GB):
///
/// * At a fixed 64-thread budget the *split* was worth 5.75x — 11.04 min with
///   every thread on mapping and no parallel decode, versus 1.92 min at
///   32 mapping / 32 decode. With one file pair the serial-gzip path caps at two
///   inflate streams no matter how many mapping threads wait on them.
/// * The knee is sharp then flat: 8 -> 16 decoder threads was worth 32%,
///   16 -> 32 worth 2%. Under-provisioning decode is expensive;
///   over-provisioning merely wastes threads.
/// * Doubling the whole budget 64 -> 128 bought only 9%.
///
/// That established the *budget* — the most decode may cost. It does not
/// establish a good starting point: measured across three workloads the knee
/// sits at 16 workers/file (Flex), 4 (bulk ERR3239276), and effectively 0
/// (PBMC on a 96 M-k-mer transcriptome, where mapping is so much costlier per
/// read that decode demand nearly vanishes). Since every curve is flat past its
/// knee, over-provisioning wastes threads without buying time. So this returns
/// a generous ceiling, a bounded total, and a deliberately small starting
/// point; `DecodeBudget` grows toward the knee only on demonstrated starvation.
///
/// `num_files` matters because throughput on the serial path scales with file
/// count (one inflate stream per file), so many inputs already supply
/// parallelism and need proportionally less help per file. Whether the parallel
/// decoder is used at all is decided by the threads-per-file ratio; see
/// [`MIN_THREADS_PER_FILE`].
pub fn plan_thread_budget(map_threads: usize, num_files: usize) -> ThreadBudget {
    plan_thread_budget_for(map_threads, num_files, ConsumerCost::Unknown)
}

/// How expensive one read is to consume, which decides how hungry the mapping
/// threads are for decoded bytes.
///
/// Callers that map the same reads two ways need this: salmon's selective
/// alignment consumes roughly a quarter of what sketch mode does per unit time
/// -- measured at 256 s of mapping-only CPU against 60 s for the same 8M reads
/// -- so the same input and thread count can be decode-bound in one mode and
/// comfortably fed in the other. Passing the mode lets them land differently
/// instead of sharing one answer that is wrong for one of them.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ConsumerCost {
    /// No hint; the threads-per-file rule alone decides. What piscem's own
    /// mapping paths use, since they have only one mode.
    #[default]
    Unknown,
    /// Cheap per read, so decode is reached sooner -- pseudoalignment/sketch.
    Light,
    /// Expensive per read, so serial decode keeps up further -- selective
    /// alignment, or mapping against a very large index.
    Heavy,
}

/// [`plan_thread_budget`] with a hint about how costly consumption is.
pub fn plan_thread_budget_for(
    map_threads: usize,
    num_files: usize,
    consumer_cost: ConsumerCost,
) -> ThreadBudget {
    if let Ok(v) = std::env::var("PISCEM_RAPIDGZIP_THREADS")
        && let Ok(per_file) = v.parse::<usize>()
    {
        // Explicit override: honour it exactly, for reproducing measurements.
        return ThreadBudget {
            decode_budget: per_file.saturating_mul(num_files.max(1)).max(1),
            per_file_ceiling: per_file.max(1),
            initial_per_file: per_file.max(1),
            parallel_gzip: true,
        };
    }

    let files = num_files.max(1);

    // Heavy consumers reach decode-bound later, so require proportionally more
    // threads per file before paying for the parallel decoder; light ones reach
    // it sooner. `Unknown` leaves the measured threshold alone.
    let threshold = match consumer_cost {
        ConsumerCost::Unknown => MIN_THREADS_PER_FILE,
        ConsumerCost::Light => MIN_THREADS_PER_FILE,
        ConsumerCost::Heavy => MIN_THREADS_PER_FILE * 2,
    };

    // Serial per-file streams already supply enough for the mapping threads.
    // This mirrors `calibrate::forced_choice`'s `AmpleSerialSupply` arm; the
    // per-file `InputKind` check lives in `open_input`, since a run can mix
    // regular files and FIFOs.
    if map_threads / files < threshold {
        return ThreadBudget {
            decode_budget: 0,
            per_file_ceiling: 0,
            initial_per_file: 0,
            parallel_gzip: false,
        };
    }
    // Half the mapping allotment, i.e. a third of the combined footprint --
    // but never more than the cores this process can actually run on.
    // `-t` is a request, not a guarantee: under a cpuset or a cgroup CPU quota
    // the process may hold far fewer cores than `-t` names, and budgeting
    // decode threads against the larger number just oversubscribes what is
    // left for mapping. (Bound suggested by salmon PR #1072.)
    let usable_cores = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(map_threads);
    let decode_budget = (map_threads / 2).clamp(1, 64).min(usable_cores.max(1));
    // Ceiling per decoder is deliberately above the fair share: a decoder that
    // is genuinely starved may exceed its slice when its peers are idle, and
    // the supervisor claws it back. Capped because past ~16 the knee is flat.
    let per_file_ceiling = (decode_budget.div_ceil(files) * 2).clamp(2, 16);

    ThreadBudget {
        decode_budget,
        per_file_ceiling,
        parallel_gzip: true,
        // Start at one worker per decoder — cost-equivalent to the serial path
        // — and let the supervisor buy growth with demonstrated starvation.
        // This removes any need to predict decode demand up front, which is
        // machine- and index-specific: measured per-mapping-thread consumption
        // ranged 0.064 GB/s (96.3 M-k-mer transcriptome) to 0.43 GB/s (1.5 M
        // probe panel), a 6.7x spread no constant could span.
        //
        // Verified not to regress the worst case: on the full Flex archive
        // (single pair, knee at 16 workers/file, where under-provisioning had
        // cost 48%), starting at 1 reached peak 16 and ran 114.29 s versus
        // 117.81 s starting at 2.
        initial_per_file: env_usize("PISCEM_DECODE_INITIAL", 1).min(per_file_ceiling).max(1),
    }
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
    #[test]
    fn test_open_concatenated_readers_empty() {
        let result = super::open_concatenated_readers(&[]);
        assert!(result.is_err());
    }
}
