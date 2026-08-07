//! Threading infrastructure — shared state for the mapping pipeline.
//!
//! Provides `OutputInfo` (mutex-guarded RAD file + chunk counter) and
//! `MappingStats` (atomic counters) used by the paraseq-based processors
//! in `crate::mapping::processors`.

use std::fs::File;
use std::io::BufWriter;
use std::sync::Mutex;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};

use crossbeam::utils::CachePadded;

// ---------------------------------------------------------------------------
// ThreadConfig
// ---------------------------------------------------------------------------

/// Threading configuration.
#[derive(Debug, Clone, Copy)]
pub struct ThreadConfig {
    /// Mapping threads, as given by `-t`.
    pub threads: usize,
    /// The user's gzip-decoder request, from `--decoder`.
    ///
    /// Lives here rather than with the mapping options because it is a
    /// statement about how to spend threads on I/O, which is what this struct
    /// is for. It reached `plan_thread_budget` via `MappingOpts` for one
    /// revision, which put a decision about decompression inside a bag of
    /// alignment parameters.
    pub decoder: crate::io::calibrate::DecoderPreference,
}

impl Default for ThreadConfig {
    fn default() -> Self {
        Self {
            threads: 16,
            decoder: crate::io::calibrate::DecoderPreference::default(),
        }
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
    /// Mutex-guarded unmapped barcode count file (scRNA/scATAC only).
    pub unmapped_bc_file: Option<Mutex<BufWriter<File>>>,
}

// ---------------------------------------------------------------------------
// MappingStats
// ---------------------------------------------------------------------------

/// Thread-safe mapping statistics.
pub struct MappingStats {
    pub num_reads: AtomicU64,
    pub num_mapped: AtomicU64,
    pub num_poisoned: AtomicU64,
    /// Time the mapping threads spent actually mapping.
    ///
    /// Lives here because it is the one struct already shared by every
    /// processor clone, so the meter needs no extra plumbing through three
    /// constructors. Read by the thread broker to decide whether the mapping
    /// threads are starved; see `crate::io::broker`.
    pub busy: std::sync::Arc<thread_broker::BusyMeter>,
    /// scATAC's current in-batch publication cadence.
    ///
    /// Other modalities use `BusyMeter`'s constant default. This indirection is
    /// read once per paraseq callback, never once per record, so the controller
    /// can make scATAC fine grained only while a decision is open without
    /// putting an atomic load in its hot loop.
    broker_progress_flush_every: std::sync::Arc<AtomicU64>,
    /// Per-processor completed-item counters used by scATAC calibration.
    /// Writers never share a cache line; the single controller reader sums the
    /// small registry only at its sampling cadence.
    broker_progress: std::sync::Arc<BrokerProgress>,
    /// Callback setup included in `busy`, and output/merge work deliberately
    /// excluded from the scalable allocation signal. Kept separately so that
    /// exclusion is measurable rather than assumed negligible.
    pub components: ConsumerComponentStats,
}

#[derive(Default)]
pub(crate) struct BrokerProgress {
    shards: Mutex<Vec<std::sync::Arc<CachePadded<AtomicU64>>>>,
}

pub(crate) struct BrokerProgressShard(std::sync::Arc<CachePadded<AtomicU64>>);

impl BrokerProgress {
    fn register(&self) -> BrokerProgressShard {
        let shard = std::sync::Arc::new(CachePadded::new(AtomicU64::new(0)));
        self.shards
            .lock()
            .unwrap()
            .push(std::sync::Arc::clone(&shard));
        BrokerProgressShard(shard)
    }

    pub(crate) fn items(&self) -> u64 {
        self.shards
            .lock()
            .unwrap()
            .iter()
            .map(|shard| shard.load(Ordering::Relaxed))
            .sum()
    }
}

impl BrokerProgressShard {
    pub(crate) fn counter(&self) -> &AtomicU64 {
        &self.0
    }
}

#[derive(Default)]
pub struct ConsumerComponentStats {
    callback_setup_nanos: AtomicU64,
    output_flush_nanos: AtomicU64,
}

impl ConsumerComponentStats {
    pub fn add_callback_setup(&self, elapsed: std::time::Duration) {
        self.callback_setup_nanos
            .fetch_add(elapsed.as_nanos() as u64, Ordering::Relaxed);
    }

    pub fn add_output_flush(&self, elapsed: std::time::Duration) {
        self.output_flush_nanos
            .fetch_add(elapsed.as_nanos() as u64, Ordering::Relaxed);
    }

    pub fn snapshot(&self) -> (u64, u64) {
        (
            self.callback_setup_nanos.load(Ordering::Relaxed),
            self.output_flush_nanos.load(Ordering::Relaxed),
        )
    }
}

impl MappingStats {
    /// Create zeroed stats.
    pub fn new() -> Self {
        Self {
            num_reads: AtomicU64::new(0),
            num_mapped: AtomicU64::new(0),
            num_poisoned: AtomicU64::new(0),
            busy: std::sync::Arc::new(thread_broker::BusyMeter::new()),
            broker_progress_flush_every: std::sync::Arc::new(AtomicU64::new(
                thread_broker::DEFAULT_FLUSH_EVERY,
            )),
            broker_progress: std::sync::Arc::new(BrokerProgress::default()),
            components: ConsumerComponentStats::default(),
        }
    }

    pub fn broker_progress_flush_every(&self) -> u64 {
        self.broker_progress_flush_every.load(Ordering::Acquire)
    }

    pub fn set_broker_progress_flush_every(&self, flush_every: u64) {
        self.broker_progress_flush_every
            .store(flush_every.max(1), Ordering::Release);
    }

    pub(crate) fn broker_progress_control(&self) -> std::sync::Arc<AtomicU64> {
        std::sync::Arc::clone(&self.broker_progress_flush_every)
    }

    pub(crate) fn broker_progress(&self) -> std::sync::Arc<BrokerProgress> {
        std::sync::Arc::clone(&self.broker_progress)
    }

    pub(crate) fn register_broker_progress_shard(&self) -> BrokerProgressShard {
        self.broker_progress.register()
    }

    /// Get summary values.
    pub fn summary(&self) -> (u64, u64, u64) {
        (
            self.num_reads.load(Ordering::Relaxed),
            self.num_mapped.load(Ordering::Relaxed),
            self.num_poisoned.load(Ordering::Relaxed),
        )
    }

    pub fn measurement_snapshot(&self) -> ConsumerMeasurement {
        let (callback_setup_nanos, output_flush_nanos) = self.components.snapshot();
        ConsumerMeasurement {
            busy_nanos: self.busy.nanos(),
            callback_setup_nanos,
            output_flush_nanos,
        }
    }

    pub fn log_measurement(
        &self,
        start: ConsumerMeasurement,
        wall: std::time::Duration,
        budget: usize,
    ) -> ConsumerMeasurement {
        let measurement = self.measurement_snapshot().delta(start);
        let capacity_nanos = (wall.as_nanos() as u64).saturating_mul(budget as u64);
        let modeled_and_excluded = measurement
            .busy_nanos
            .saturating_add(measurement.output_flush_nanos);
        tracing::info!(
            modeled_busy_seconds = measurement.busy_nanos as f64 / 1e9,
            callback_setup_seconds = measurement.callback_setup_nanos as f64 / 1e9,
            excluded_output_seconds = measurement.output_flush_nanos as f64 / 1e9,
            wall_seconds = wall.as_secs_f64(),
            budget,
            modeled_utilization = if capacity_nanos > 0 {
                measurement.busy_nanos as f64 / capacity_nanos as f64
            } else {
                0.0
            },
            excluded_output_fraction = if modeled_and_excluded > 0 {
                measurement.output_flush_nanos as f64 / modeled_and_excluded as f64
            } else {
                0.0
            },
            "consumer measurement components"
        );
        measurement
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, serde::Serialize)]
pub struct ConsumerMeasurement {
    pub busy_nanos: u64,
    pub callback_setup_nanos: u64,
    pub output_flush_nanos: u64,
}

impl ConsumerMeasurement {
    pub fn delta(self, previous: Self) -> Self {
        Self {
            busy_nanos: self.busy_nanos.saturating_sub(previous.busy_nanos),
            callback_setup_nanos: self
                .callback_setup_nanos
                .saturating_sub(previous.callback_setup_nanos),
            output_flush_nanos: self
                .output_flush_nanos
                .saturating_sub(previous.output_flush_nanos),
        }
    }
}

impl Default for MappingStats {
    fn default() -> Self {
        Self::new()
    }
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

    #[test]
    fn broker_progress_is_sharded_and_cumulative() {
        let stats = MappingStats::new();
        let left = stats.register_broker_progress_shard();
        let right = stats.register_broker_progress_shard();
        left.counter().fetch_add(64, Ordering::Relaxed);
        right.counter().fetch_add(32, Ordering::Relaxed);
        assert_eq!(stats.broker_progress.items(), 96);
    }
}
