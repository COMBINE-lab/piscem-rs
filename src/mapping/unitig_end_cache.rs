//! Concurrent unitig-end cache for avoiding redundant dictionary lookups.
//!
//! When a streaming query reaches the end of a unitig, the next k-mer
//! requires a full (expensive) dictionary search. If another thread has
//! already searched for the same canonical k-mer, we can reuse the cached
//! result. This is especially effective for shared unitig junctions.
//!
//! Port of C++ `boost::concurrent_flat_map<uint64_t, lookup_result>` from
//! `piscem-cpp/include/streaming_query.hpp`.

use dashmap::DashMap;
use sshash_lib::LookupResult;

/// Cached lookup result with orientation metadata.
#[derive(Clone, Debug)]
struct CachedLookup {
    result: LookupResult,
    /// Whether the forward k-mer was canonical at insertion time.
    fw_was_canonical: bool,
}

/// Thread-safe concurrent cache for unitig-end k-mer lookups.
///
/// Shared across all worker threads via `&UnitigEndCache`. Insertions
/// are bounded by `capacity` — once the cache is full, no new entries
/// are inserted (existing entries remain accessible).
pub struct UnitigEndCache {
    map: DashMap<u64, CachedLookup>,
    capacity: usize,
}

impl UnitigEndCache {
    /// Create a new cache with the given maximum capacity.
    pub fn new(capacity: usize) -> Self {
        Self {
            map: DashMap::with_capacity(capacity.min(1 << 20)),
            capacity,
        }
    }

    /// Try to retrieve a cached lookup for the given canonical k-mer hash.
    ///
    /// If the stored orientation doesn't match the current query's
    /// `fw_is_canonical`, the result's orientation is flipped.
    pub fn get(&self, canonical_hash: u64, fw_is_canonical: bool) -> Option<LookupResult> {
        self.map.get(&canonical_hash).map(|entry| {
            let cached = entry.value();
            let mut result = cached.result.clone();
            if fw_is_canonical != cached.fw_was_canonical {
                // Flip orientation: forward ↔ reverse complement
                result.kmer_orientation = -result.kmer_orientation;
            }
            result
        })
    }

    /// Insert a lookup result into the cache.
    ///
    /// Does nothing if the cache is at capacity.
    pub fn insert(&self, canonical_hash: u64, result: &LookupResult, fw_is_canonical: bool) {
        if self.map.len() >= self.capacity {
            return;
        }
        self.map.insert(
            canonical_hash,
            CachedLookup {
                result: result.clone(),
                fw_was_canonical: fw_is_canonical,
            },
        );
    }

    /// Number of entries currently in the cache.
    pub fn len(&self) -> usize {
        self.map.len()
    }

    /// Whether the cache is empty.
    pub fn is_empty(&self) -> bool {
        self.map.is_empty()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn make_result(string_id: u64, kmer_id: u64, orientation: i8) -> LookupResult {
        LookupResult {
            kmer_id,
            kmer_id_in_string: 5,
            kmer_offset: 0,
            kmer_orientation: orientation,
            string_id,
            string_begin: 0,
            string_end: 100,
            minimizer_found: true,
        }
    }

    #[test]
    fn test_cache_insert_and_retrieve() {
        let cache = UnitigEndCache::new(100);

        let result = make_result(42, 1000, 1); // forward
        cache.insert(0xDEAD, &result, true); // fw_is_canonical=true

        // Same orientation → no flip
        let got = cache.get(0xDEAD, true).unwrap();
        assert_eq!(got.kmer_orientation, 1);
        assert_eq!(got.string_id, 42);
    }

    #[test]
    fn test_cache_orientation_flip() {
        let cache = UnitigEndCache::new(100);

        let result = make_result(42, 1000, 1); // forward
        cache.insert(0xBEEF, &result, true); // stored with fw_is_canonical=true

        // Different fw_is_canonical → orientation flipped
        let got = cache.get(0xBEEF, false).unwrap();
        assert_eq!(got.kmer_orientation, -1); // flipped from 1 to -1
    }

    #[test]
    fn test_cache_capacity_limit() {
        let cache = UnitigEndCache::new(3);

        for i in 0..5u64 {
            let result = make_result(i, i * 10, 1);
            cache.insert(i, &result, true);
        }

        // Only first 3 should be present
        assert!(cache.len() <= 3);
        assert!(cache.get(0, true).is_some());
        assert!(cache.get(1, true).is_some());
        assert!(cache.get(2, true).is_some());
    }

    #[test]
    fn test_cache_miss_falls_through() {
        let cache = UnitigEndCache::new(100);
        assert!(cache.get(0xDEAD, true).is_none());
        assert!(cache.is_empty());
    }
}
