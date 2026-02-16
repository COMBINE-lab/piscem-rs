//! Threading infrastructure — shared state for the mapping pipeline.
//!
//! Provides `OutputInfo` (mutex-guarded RAD file + chunk counter) and
//! `MappingStats` (atomic counters) used by the paraseq-based processors
//! in `crate::mapping::processors`.

use std::fs::File;
use std::io::BufWriter;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::Mutex;

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
