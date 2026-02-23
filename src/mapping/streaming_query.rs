//! Thin wrapper around the sshash-rs streaming query engine.
//!
//! `PiscemStreamingQuery` adapts the sshash-rs `StreamingQueryEngine` for use
//! in the piscem mapping pipeline. It provides the same `lookup()` / `reset()`
//! interface, tracks search statistics, and integrates the unitig-end cache.
//!
//! Corresponds to the C++ `piscem::streaming_query`.

use sshash_lib::streaming_query::StreamingQueryEngine;
use sshash_lib::{Dictionary, Kmer, KmerBits, LookupResult};

use crate::mapping::kmer_value::CanonicalKmer;
use crate::mapping::unitig_end_cache::UnitigEndCache;

/// Streaming k-mer query engine for the piscem mapping pipeline.
///
/// Wraps sshash-rs `StreamingQueryEngine`, which performs efficient
/// k-mer lookups by extending along unitigs when possible (avoiding
/// full dictionary searches for consecutive k-mers on the same unitig).
///
/// Optionally integrates with a shared `UnitigEndCache` to avoid
/// redundant lookups at unitig boundaries.
pub struct PiscemStreamingQuery<'a, const K: usize>
where
    Kmer<K>: KmerBits,
{
    engine: StreamingQueryEngine<'a, K>,
    k: usize,
    /// Set to `true` when we reach a unitig boundary, triggering
    /// cache lookup/insert on the next query.
    cache_end: bool,
    /// Optional shared unitig-end cache.
    cache: Option<&'a UnitigEndCache>,
    /// Number of cache hits (for statistics).
    num_cache_hits: u64,
    /// Position of the previous query in the read (for auto-reset detection).
    /// When the next lookup is not at `prev_query_pos + 1`, the engine is
    /// automatically reset, forcing a full k-mer parse. This mirrors the
    /// C++ `m_prev_query_offset` tracking and allows the sshash-rs engine
    /// to use its fast incremental k-mer update path for consecutive lookups.
    prev_query_pos: i32,
}

impl<'a, const K: usize> PiscemStreamingQuery<'a, K>
where
    Kmer<K>: KmerBits,
{
    /// Create a new streaming query engine from a dictionary.
    pub fn new(dict: &'a Dictionary) -> Self {
        Self {
            k: dict.k(),
            engine: dict.create_streaming_query::<K>(),
            cache_end: false,
            cache: None,
            num_cache_hits: 0,
            prev_query_pos: i32::MIN,
        }
    }

    /// Create a new streaming query engine with a shared unitig-end cache.
    pub fn with_cache(dict: &'a Dictionary, cache: &'a UnitigEndCache) -> Self {
        Self {
            k: dict.k(),
            engine: dict.create_streaming_query::<K>(),
            cache_end: false,
            cache: Some(cache),
            num_cache_hits: 0,
            prev_query_pos: i32::MIN,
        }
    }

    /// Reset the query state. Call this between read sequences so that the
    /// engine doesn't try to extend from a previous sequence's unitig.
    #[inline]
    pub fn reset(&mut self) {
        self.engine.reset();
        self.cache_end = false;
        self.prev_query_pos = i32::MIN;
    }

    /// Look up a k-mer at a specific read position.
    ///
    /// Tracks position to auto-detect non-consecutive lookups. If `read_pos`
    /// is exactly `prev_query_pos + 1`, the sshash-rs engine reuses its
    /// incremental k-mer state (fast). Otherwise, the engine is reset and
    /// the k-mer is parsed from scratch.
    ///
    /// This mirrors the C++ `m_prev_query_offset` tracking in
    /// `piscem::streaming_query`.
    ///
    /// If the unitig-end cache is enabled and we're at a unitig boundary,
    /// checks the cache first. On a cache miss, performs a full lookup and
    /// caches the result for future threads.
    #[inline]
    pub fn lookup_at(&mut self, kmer_bytes: &[u8], read_pos: i32) -> LookupResult {
        // Auto-detect non-consecutive position and reset engine if needed.
        // This recovers the incremental k-mer update path for consecutive
        // lookups (stride=1) while correctly handling jumps.
        if read_pos != self.prev_query_pos.wrapping_add(1) {
            self.engine.reset();
            self.cache_end = false;
        }
        self.prev_query_pos = read_pos;

        // Try cache if we're at a unitig boundary
        if self.cache_end
            && let Some(cache) = self.cache
        {
            let kmer = Kmer::<K>::from_ascii_unchecked(kmer_bytes);
            let canonical = kmer.canonical();
            let canonical_hash =
                CanonicalKmer::new(<Kmer<K> as KmerBits>::to_u64(canonical.bits()));
            let fw_is_canonical = kmer.bits() == canonical.bits();

            if let Some(result) = cache.get(canonical_hash, fw_is_canonical) {
                self.num_cache_hits += 1;
                self.cache_end = false;
                // Reset engine state since we bypassed it
                self.engine.reset();
                return result;
            }
        }

        let was_cache_end = self.cache_end;
        self.cache_end = false;

        // Full lookup via streaming query engine
        let result = self.engine.lookup(kmer_bytes);

        // If the lookup succeeded and we were at a unitig boundary, cache it
        if was_cache_end && result.is_found()
            && let Some(cache) = self.cache
        {
            let kmer = Kmer::<K>::from_ascii_unchecked(kmer_bytes);
            let canonical = kmer.canonical();
            let canonical_hash =
                CanonicalKmer::new(<Kmer<K> as KmerBits>::to_u64(canonical.bits()));
            let fw_is_canonical = kmer.bits() == canonical.bits();
            cache.insert(canonical_hash, &result, fw_is_canonical);
        }

        // Track remaining bases to detect unitig boundaries
        if result.is_found() {
            let direction = result.kmer_orientation;
            let remaining = if direction > 0 {
                // Forward: remaining = string_length - (kmer_id_in_string + k)
                let slen = result.string_length();
                slen.saturating_sub(result.kmer_id_in_string + self.k as u64)
            } else {
                // Backward: remaining = kmer_id_in_string
                result.kmer_id_in_string
            };
            if remaining == 0 {
                self.cache_end = true;
            }
        }

        result
    }

    /// Look up a k-mer without position tracking (always resets engine).
    ///
    /// This is a convenience method for callers that don't track read position.
    /// For optimal performance in the mapping pipeline, use `lookup_at()` instead.
    #[inline]
    pub fn lookup(&mut self, kmer_bytes: &[u8]) -> LookupResult {
        self.engine.reset();
        self.prev_query_pos = i32::MIN;
        self.lookup_at(kmer_bytes, 0)
    }

    /// Number of full dictionary searches performed (expensive path).
    #[inline]
    pub fn num_searches(&self) -> u64 {
        self.engine.num_searches()
    }

    /// Number of unitig extensions performed (cheap path).
    #[inline]
    pub fn num_extensions(&self) -> u64 {
        self.engine.num_extensions()
    }

    /// Number of cache hits.
    #[inline]
    pub fn num_cache_hits(&self) -> u64 {
        self.num_cache_hits
    }
}
