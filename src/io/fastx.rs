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
/// Swept and confirmed on the Flex pair at `-t 32` through `io::probe`, which
/// runs the real pipeline, at three consumer costs. 4096 and 16384 form a
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
    // a cheap starting point the supervisor could ratchet up. It was instead a
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

    // `ceiling == 0` selects the serial path explicitly (see
    // the probe's verdict), independent of whether the feature is on.
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

/// How a run intends to spend its thread budget.
#[derive(Debug, Clone, Copy)]
pub struct ThreadBudget {
    /// Threads that will map. Together with `decode_budget` this sums to the
    /// user's `-t`.
    pub map_threads: usize,
    /// Total decoder workers all files may hold at once.
    pub decode_budget: usize,
    /// Immutable per-decoder ceiling. Generous on purpose: it is a lazy cap,
    /// and `DecodeBudget` enforces the real limit at runtime.
    pub per_file_ceiling: usize,
    /// Whether to use the parallel gzip decoder at all. When false, gzip input
    /// takes the serial niffler path and no decoder threads are spawned.
    pub parallel_gzip: bool,
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
pub fn plan_thread_budget(total_threads: usize, num_files: usize) -> ThreadBudget {
    let files = num_files.max(1);
    let total = total_threads.max(1);

    if let Ok(v) = std::env::var("PISCEM_RAPIDGZIP_THREADS")
        && let Ok(per_file) = v.parse::<usize>()
    {
        // Explicit override: honour it exactly, for reproducing measurements.
        let decode = per_file.saturating_mul(files).max(1);
        return ThreadBudget {
            map_threads: total.saturating_sub(decode).max(1),
            decode_budget: decode,
            per_file_ceiling: per_file.max(1),
            parallel_gzip: true,
        };
    }

    // `-t` is a request, not a guarantee: under a cpuset or a cgroup CPU quota
    // the process may hold far fewer cores than it names, and budgeting against
    // the larger number just oversubscribes. (Bound suggested by salmon #1072.)
    let usable = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(total);
    let total = total.min(usable.max(1));

    let decode_budget = (total / 8).clamp(1, total.saturating_sub(1).max(1));
    let map_threads = total.saturating_sub(decode_budget).max(1);

    // The per-file ceiling is now the *pool's* immutable maximum, so it is the
    // whole budget rather than a share of it. `DecoderPool::set_worker_limit` is
    // refused above this, and a refusal silently desyncs the broker's accounting
    // from the pool's -- so the ceiling must never be the thing that binds.
    ThreadBudget {
        map_threads,
        decode_budget,
        per_file_ceiling: total,
        parallel_gzip: true,
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
