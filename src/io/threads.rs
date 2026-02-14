//! Threading infrastructure — producer-consumer pipeline.
//!
//! Provides the mapping pipeline orchestrator that reads FASTX chunks,
//! dispatches them to worker threads via a bounded crossbeam channel,
//! and collects output.
//!
//! Uses crossbeam scoped threads for natural lifetime management: worker
//! threads can borrow the index and other shared state without `Arc`.

use std::fs::File;
use std::io::BufWriter;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::Mutex;

use anyhow::Result;
use crossbeam::channel;

use crate::io::fastx::{FastxSource, FastxTripleSource, ReadChunk, ReadTripletChunk};

// ---------------------------------------------------------------------------
// ThreadConfig
// ---------------------------------------------------------------------------

/// Threading configuration.
#[derive(Debug, Clone, Copy)]
pub struct ThreadConfig {
    pub threads: usize,
}

impl Default for ThreadConfig {
    fn default() -> Self {
        Self { threads: 16 }
    }
}

// ---------------------------------------------------------------------------
// OutputInfo
// ---------------------------------------------------------------------------

/// Shared output state for the mapping pipeline.
pub struct OutputInfo {
    /// Number of chunks processed so far.
    pub num_chunks: AtomicUsize,
    /// Mutex-guarded RAD output file.
    pub rad_file: Mutex<BufWriter<File>>,
}

// ---------------------------------------------------------------------------
// MappingStats
// ---------------------------------------------------------------------------

/// Thread-safe mapping statistics.
pub struct MappingStats {
    pub num_reads: AtomicU64,
    pub num_mapped: AtomicU64,
    pub num_poisoned: AtomicU64,
}

impl MappingStats {
    /// Create zeroed stats.
    pub fn new() -> Self {
        Self {
            num_reads: AtomicU64::new(0),
            num_mapped: AtomicU64::new(0),
            num_poisoned: AtomicU64::new(0),
        }
    }

    /// Get summary values.
    pub fn summary(&self) -> (u64, u64, u64) {
        (
            self.num_reads.load(Ordering::Relaxed),
            self.num_mapped.load(Ordering::Relaxed),
            self.num_poisoned.load(Ordering::Relaxed),
        )
    }
}

impl Default for MappingStats {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// run_mapping_pipeline
// ---------------------------------------------------------------------------

/// Run the mapping pipeline.
///
/// - Main thread reads FASTX chunks, sends through a bounded channel
/// - N workers receive chunks, process reads, flush RAD output under mutex
/// - Single-thread mode: bypass channel, process inline
///
/// The `worker_fn` closure is where protocol-specific logic goes. Each
/// protocol's CLI command constructs the closure that creates per-thread
/// `MappingCache`, `HitSearcher`, `PiscemStreamingQuery`, `PoisonState`,
/// processes reads, and flushes RAD output.
pub fn run_mapping_pipeline<F>(
    mut fastx: FastxSource,
    config: ThreadConfig,
    output: &OutputInfo,
    stats: &MappingStats,
    worker_fn: F,
) -> Result<()>
where
    F: Fn(ReadChunk, &OutputInfo, &MappingStats) + Send + Sync,
{
    let num_threads = config.threads.max(1);

    // Always use producer-consumer pattern: 1 dedicated producer thread for
    // FASTQ parsing/decompression + N worker threads for mapping.
    // This matches C++ piscem which always uses a separate parser thread.
    // At t=1: 1 producer + 1 worker = 2 OS threads (parallel I/O + mapping).
    let (sender, receiver) = channel::bounded::<ReadChunk>(num_threads * 2);

    let worker_ref = &worker_fn;
    crossbeam::scope(|scope| {
        // Spawn worker threads
        for _ in 0..num_threads {
            let recv = receiver.clone();
            scope.spawn(move |_| {
                while let Ok(chunk) = recv.recv() {
                    worker_ref(chunk, output, stats);
                }
            });
        }
        // Drop the extra receiver clone so workers will exit when sender is dropped.
        drop(receiver);

        // Producer thread: reads FASTX chunks and sends to workers.
        scope.spawn(move |_| {
            let mut chunk = Vec::new();
            loop {
                match fastx.next_chunk(&mut chunk) {
                    Ok(true) => {
                        let batch = std::mem::take(&mut chunk);
                        if sender.send(batch).is_err() {
                            break; // Workers have shut down
                        }
                    }
                    Ok(false) => break, // EOF
                    Err(e) => {
                        tracing::error!("Error reading FASTX: {}", e);
                        break;
                    }
                }
            }
            // sender dropped here, signaling workers to exit.
        });
    })
    .map_err(|e| anyhow::anyhow!("thread panicked: {:?}", e))?;

    Ok(())
}

// ---------------------------------------------------------------------------
// run_mapping_pipeline_triple (for scATAC triple-file input)
// ---------------------------------------------------------------------------

/// Run the mapping pipeline with triple-file FASTQ input (R1 + barcode + R2).
///
/// Same producer-consumer pattern as `run_mapping_pipeline`, but uses
/// `FastxTripleSource` and `ReadTripletChunk`.
pub fn run_mapping_pipeline_triple<F>(
    mut fastx: FastxTripleSource,
    config: ThreadConfig,
    output: &OutputInfo,
    stats: &MappingStats,
    worker_fn: F,
) -> Result<()>
where
    F: Fn(ReadTripletChunk, &OutputInfo, &MappingStats) + Send + Sync,
{
    let num_threads = config.threads.max(1);
    let (sender, receiver) = channel::bounded::<ReadTripletChunk>(num_threads * 2);

    let worker_ref = &worker_fn;
    crossbeam::scope(|scope| {
        // Spawn worker threads
        for _ in 0..num_threads {
            let recv = receiver.clone();
            scope.spawn(move |_| {
                while let Ok(chunk) = recv.recv() {
                    worker_ref(chunk, output, stats);
                }
            });
        }
        drop(receiver);

        // Producer thread
        scope.spawn(move |_| {
            let mut chunk = Vec::new();
            loop {
                match fastx.next_chunk(&mut chunk) {
                    Ok(true) => {
                        let batch = std::mem::take(&mut chunk);
                        if sender.send(batch).is_err() {
                            break;
                        }
                    }
                    Ok(false) => break,
                    Err(e) => {
                        tracing::error!("Error reading FASTX: {}", e);
                        break;
                    }
                }
            }
        });
    })
    .map_err(|e| anyhow::anyhow!("thread panicked: {:?}", e))?;

    Ok(())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_thread_config_default() {
        let config = ThreadConfig::default();
        assert_eq!(config.threads, 16);
    }

    #[test]
    fn test_mapping_stats_new() {
        let stats = MappingStats::new();
        let (r, m, p) = stats.summary();
        assert_eq!(r, 0);
        assert_eq!(m, 0);
        assert_eq!(p, 0);
    }

    #[test]
    fn test_mapping_stats_atomic_ops() {
        let stats = MappingStats::new();
        stats.num_reads.fetch_add(100, Ordering::Relaxed);
        stats.num_mapped.fetch_add(80, Ordering::Relaxed);
        stats.num_poisoned.fetch_add(5, Ordering::Relaxed);
        let (r, m, p) = stats.summary();
        assert_eq!(r, 100);
        assert_eq!(m, 80);
        assert_eq!(p, 5);
    }
}
