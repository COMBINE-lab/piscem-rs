//! Mapping cache — per-read scratch state for the mapping kernel.
//!
//! Port of C++ `mapping_cache_info<T>`. Generic over the hit info type
//! (`SketchHitInfoSimple` or `SketchHitInfoChained`).
//!
//! Note: The streaming query and hit searcher are NOT owned here — they
//! borrow the index and are created as separate per-thread locals alongside
//! the cache.

use std::collections::{HashMap, HashSet};

use nohash_hasher::BuildNoHashHasher;

use crate::mapping::hits::{MappingType, SimpleHit, SketchHitInfo};

// ---------------------------------------------------------------------------
// MappingCache
// ---------------------------------------------------------------------------

/// Per-read mapping state, generic over the hit accumulator type.
///
/// Corresponds to C++ `mapping_cache_info<sketch_hit_info_t>`.
pub struct MappingCache<S: SketchHitInfo> {
    /// How the read mapped.
    pub map_type: MappingType,
    /// Map from reference target ID → per-target hit accumulator.
    pub hit_map: HashMap<u32, S, BuildNoHashHasher<u32>>,
    /// Final accepted hit list (after filtering).
    pub accepted_hits: Vec<SimpleHit>,
    /// Maximum number of reference occurrences before a k-mer is considered
    /// too ambiguous (first pass threshold).
    pub max_hit_occ: usize,
    /// Recovery threshold: retry with `min_occ` if all k-mers exceeded
    /// `max_hit_occ` but `min_occ < max_hit_occ_recover`.
    pub max_hit_occ_recover: usize,
    /// Whether occurrence recovery is enabled.
    pub attempt_occ_recover: bool,
    /// Maximum number of accepted mappings before discarding the read.
    pub max_read_occ: usize,
    /// K-mer size.
    pub k: usize,
    /// Maximum equivalence class cardinality for ambiguous hit filtering.
    pub max_ec_card: u32,
    /// Whether any k-mers produced index matches (even if too ambiguous).
    pub has_matching_kmers: bool,
    /// Indices of raw hits that exceeded the occurrence threshold (for EC filtering).
    pub ambiguous_hit_indices: Vec<u32>,
    /// Reusable set for tracking observed ECs during ambiguous hit filtering.
    /// Kept here to avoid per-read allocation.
    pub observed_ecs: HashSet<u64>,
}

impl<S: SketchHitInfo> MappingCache<S> {
    /// Create a new mapping cache with the given k-mer size and defaults.
    pub fn new(k: usize) -> Self {
        let max_hit_occ = 256;
        let max_hit_occ_recover = 1024;
        Self {
            map_type: MappingType::Unmapped,
            hit_map: HashMap::with_hasher(BuildNoHashHasher::default()),
            accepted_hits: Vec::new(),
            max_hit_occ,
            max_hit_occ_recover,
            attempt_occ_recover: max_hit_occ_recover > max_hit_occ,
            max_read_occ: 2500,
            k,
            max_ec_card: 4096,
            has_matching_kmers: false,
            ambiguous_hit_indices: Vec::new(),
            observed_ecs: HashSet::new(),
        }
    }

    /// Reset state between reads.
    pub fn clear(&mut self) {
        self.map_type = MappingType::Unmapped;
        self.hit_map.clear();
        self.accepted_hits.clear();
        self.has_matching_kmers = false;
        self.ambiguous_hit_indices.clear();
        self.observed_ecs.clear();
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mapping::sketch_hit_simple::SketchHitInfoSimple;

    #[test]
    fn test_cache_defaults() {
        let cache = MappingCache::<SketchHitInfoSimple>::new(31);
        assert_eq!(cache.map_type, MappingType::Unmapped);
        assert!(cache.hit_map.is_empty());
        assert!(cache.accepted_hits.is_empty());
        assert_eq!(cache.max_hit_occ, 256);
        assert_eq!(cache.max_hit_occ_recover, 1024);
        assert!(cache.attempt_occ_recover);
        assert_eq!(cache.max_read_occ, 2500);
        assert_eq!(cache.k, 31);
        assert_eq!(cache.max_ec_card, 4096);
        assert!(!cache.has_matching_kmers);
    }

    #[test]
    fn test_cache_clear() {
        let mut cache = MappingCache::<SketchHitInfoSimple>::new(31);
        cache.map_type = MappingType::SingleMapped;
        cache.has_matching_kmers = true;
        cache.accepted_hits.push(SimpleHit::default());
        cache.ambiguous_hit_indices.push(42);
        cache
            .hit_map
            .insert(0, SketchHitInfoSimple::default());

        cache.clear();

        assert_eq!(cache.map_type, MappingType::Unmapped);
        assert!(cache.hit_map.is_empty());
        assert!(cache.accepted_hits.is_empty());
        assert!(!cache.has_matching_kmers);
        assert!(cache.ambiguous_hit_indices.is_empty());
    }
}
