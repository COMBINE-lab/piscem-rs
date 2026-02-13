//! Thin wrapper around the sshash-rs streaming query engine.
//!
//! `PiscemStreamingQuery` adapts the sshash-rs `StreamingQueryEngine` for use
//! in the piscem mapping pipeline. It provides the same `lookup()` / `reset()`
//! interface, tracks search statistics, and integrates the unitig-end cache.
//!
//! Corresponds to the C++ `piscem::streaming_query`.

use sshash_lib::streaming_query::StreamingQueryEngine;
use sshash_lib::{Dictionary, Kmer, KmerBits, LookupResult};

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
        }
    }

    /// Reset the query state. Call this between read sequences so that the
    /// engine doesn't try to extend from a previous sequence's unitig.
    #[inline]
    pub fn reset(&mut self) {
        self.engine.reset();
        self.cache_end = false;
    }

    /// Look up a k-mer string in the dictionary.
    ///
    /// If the unitig-end cache is enabled and we're at a unitig boundary,
    /// checks the cache first. On a cache miss, performs a full lookup and
    /// caches the result for future threads.
    ///
    /// Returns a `LookupResult` which can be checked with `is_found()` and
    /// then passed to `ReferenceIndex::resolve_lookup()`.
    #[inline]
    pub fn lookup(&mut self, kmer_str: &str) -> LookupResult {
        // Try cache if we're at a unitig boundary
        if self.cache_end {
            if let Some(cache) = self.cache {
                if let Ok(kmer) = Kmer::<K>::from_str(kmer_str) {
                    let canonical = kmer.canonical();
                    let canonical_hash =
                        <Kmer<K> as KmerBits>::to_u64(canonical.bits());
                    let fw_is_canonical = kmer.bits() == canonical.bits();

                    if let Some(result) = cache.get(canonical_hash, fw_is_canonical) {
                        self.num_cache_hits += 1;
                        self.cache_end = false;
                        // Reset engine state since we bypassed it
                        self.engine.reset();
                        return result;
                    }
                }
            }
        }

        let was_cache_end = self.cache_end;
        self.cache_end = false;

        // Full lookup via streaming query engine
        let result = self.engine.lookup(kmer_str);

        // If the lookup succeeded and we were at a unitig boundary, cache it
        if was_cache_end && result.is_found() {
            if let Some(cache) = self.cache {
                if let Ok(kmer) = Kmer::<K>::from_str(kmer_str) {
                    let canonical = kmer.canonical();
                    let canonical_hash =
                        <Kmer<K> as KmerBits>::to_u64(canonical.bits());
                    let fw_is_canonical = kmer.bits() == canonical.bits();
                    cache.insert(canonical_hash, &result, fw_is_canonical);
                }
            }
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
