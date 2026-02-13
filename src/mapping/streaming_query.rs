//! Thin wrapper around the sshash-rs streaming query engine.
//!
//! `PiscemStreamingQuery` adapts the sshash-rs `StreamingQueryEngine` for use
//! in the piscem mapping pipeline. It provides the same `lookup()` / `reset()`
//! interface and tracks search statistics.
//!
//! Corresponds to the C++ `piscem::streaming_query` (without the unitig-end
//! cache, which is an optimization deferred to Phase 5).

use sshash_lib::streaming_query::StreamingQueryEngine;
use sshash_lib::{Dictionary, Kmer, KmerBits, LookupResult};

/// Streaming k-mer query engine for the piscem mapping pipeline.
///
/// Wraps sshash-rs `StreamingQueryEngine`, which performs efficient
/// k-mer lookups by extending along unitigs when possible (avoiding
/// full dictionary searches for consecutive k-mers on the same unitig).
pub struct PiscemStreamingQuery<'a, const K: usize>
where
    Kmer<K>: KmerBits,
{
    engine: StreamingQueryEngine<'a, K>,
}

impl<'a, const K: usize> PiscemStreamingQuery<'a, K>
where
    Kmer<K>: KmerBits,
{
    /// Create a new streaming query engine from a dictionary.
    pub fn new(dict: &'a Dictionary) -> Self {
        Self {
            engine: dict.create_streaming_query::<K>(),
        }
    }

    /// Reset the query state. Call this between read sequences so that the
    /// engine doesn't try to extend from a previous sequence's unitig.
    #[inline]
    pub fn reset(&mut self) {
        self.engine.reset();
    }

    /// Look up a k-mer string in the dictionary.
    ///
    /// Returns a `LookupResult` which can be checked with `is_found()` and
    /// then passed to `ReferenceIndex::resolve_lookup()`.
    #[inline]
    pub fn lookup(&mut self, kmer_str: &str) -> LookupResult {
        self.engine.lookup(kmer_str)
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
}
