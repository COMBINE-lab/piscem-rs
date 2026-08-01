//! Mapping cache — per-read scratch state for the mapping kernel.
//!
//! Port of C++ `mapping_cache_info<T>`. Generic over the hit info type
//! (`SketchHitInfoSimple` or `SketchHitInfoChained`).
//!
//! Note: The streaming query and hit searcher are NOT owned here — they
//! borrow the index and are created as separate per-thread locals alongside
//! the cache.

use crate::hash::{fixed_map_with_capacity, fixed_set, FixedMap, FixedSet};

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
    /// Uses AHashMap (fast hash) matching C++ ankerl::unordered_dense::map performance.
    pub hit_map: FixedMap<u32, S>,
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
    /// Lower bound of the accepted paired fragment-length window in
    /// `merge_se_mappings` (the upper bound is fixed at 2000). Default `-32`
    /// (piscem's small-dovetail tolerance); set to a more negative value
    /// (e.g. `-read_len`) to admit dovetailed fragments (salmon `--allowDovetail`).
    pub min_frag_len: i32,
    /// Orphan policy for `merge_se_mappings`. Default `false` = piscem's strict
    /// rule (orphan a pair only when the other mate had *no matching k-mers*).
    /// `true` orphans whenever exactly one mate has an accepted target,
    /// regardless of the other mate's k-mers (used by salmon's sketch path,
    /// where it improves agreement with selective alignment).
    pub relaxed_orphans: bool,
    /// K-mer size.
    pub k: usize,
    /// Maximum equivalence class cardinality for ambiguous hit filtering.
    pub max_ec_card: u32,
    /// Whether any k-mers produced index matches (even if too ambiguous).
    pub has_matching_kmers: bool,
    /// Indices of raw hits that exceeded the occurrence threshold (for EC filtering).
    pub ambiguous_hit_indices: Vec<u32>,
    /// Reusable set for tracking observed ECs during ambiguous hit filtering.
    /// Uses AHashSet (fast hash) instead of standard HashSet (SipHash) to match
    /// C++ performance with ankerl::unordered_dense::set.
    pub observed_ecs: FixedSet<u64>,
}

impl<S: SketchHitInfo> MappingCache<S> {
    /// Create a new mapping cache with the given k-mer size and defaults.
    pub fn new(k: usize) -> Self {
        let max_hit_occ = 256;
        let max_hit_occ_recover = 1024;
        // Pre-reserve capacity matching C++ (utils.hpp line 932, 966)
        let hit_map = fixed_map_with_capacity(max_hit_occ);
        Self {
            map_type: MappingType::Unmapped,
            hit_map,
            accepted_hits: Vec::with_capacity(64),
            max_hit_occ,
            max_hit_occ_recover,
            attempt_occ_recover: max_hit_occ_recover > max_hit_occ,
            max_read_occ: 2500,
            min_frag_len: -32,
            relaxed_orphans: false,
            k,
            max_ec_card: 4096,
            has_matching_kmers: false,
            ambiguous_hit_indices: Vec::new(),
            observed_ecs: fixed_set(),
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
        cache.hit_map.insert(0, SketchHitInfoSimple::default());

        cache.clear();

        assert_eq!(cache.map_type, MappingType::Unmapped);
        assert!(cache.hit_map.is_empty());
        assert!(cache.accepted_hits.is_empty());
        assert!(!cache.has_matching_kmers);
        assert!(cache.ambiguous_hit_indices.is_empty());
    }
}
