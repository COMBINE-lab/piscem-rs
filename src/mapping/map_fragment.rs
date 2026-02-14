//! Per-read mapping helpers — `map_se_fragment` / `map_pe_fragment`.
//!
//! Factors out common per-read mapping flow shared by all protocol
//! implementations.
//!
//! Port of C++ `map_se_fragment` / `map_pe_fragment` from
//! `piscem-cpp/src/pesc_sc.cpp`.

use sshash_lib::{Kmer, KmerBits};

use crate::index::reference_index::ReferenceIndex;
use crate::mapping::binning::BinPos;
use crate::mapping::cache::MappingCache;
use crate::mapping::engine::{map_read, map_read_atac};
use crate::mapping::filters::PoisonState;
use crate::mapping::hit_searcher::{HitSearcher, SkippingStrategy};
use crate::mapping::hits::FragmentEnd;
use crate::mapping::merge_pairs::{merge_se_mappings, merge_se_mappings_binned};
use crate::mapping::sketch_hit_simple::SketchHitInfoSimple;
use crate::mapping::streaming_query::PiscemStreamingQuery;

/// Map a single-end read fragment.
///
/// Returns `true` for early stop.
pub fn map_se_fragment<const K: usize>(
    seq: &[u8],
    hs: &mut HitSearcher<'_>,
    query: &mut PiscemStreamingQuery<'_, K>,
    cache_out: &mut MappingCache<SketchHitInfoSimple>,
    index: &ReferenceIndex,
    poison_state: &mut PoisonState<'_>,
    strat: SkippingStrategy,
) -> bool
where
    Kmer<K>: KmerBits,
{
    poison_state.clear();
    map_read::<K, SketchHitInfoSimple>(seq, cache_out, hs, query, index, poison_state, strat)
}

/// Map a paired-end read fragment.
///
/// Maps both ends independently, then merges results.
/// Returns `true` for early stop.
#[allow(clippy::too_many_arguments)]
pub fn map_pe_fragment<const K: usize>(
    seq1: &[u8],
    seq2: &[u8],
    hs: &mut HitSearcher<'_>,
    query: &mut PiscemStreamingQuery<'_, K>,
    cache_left: &mut MappingCache<SketchHitInfoSimple>,
    cache_right: &mut MappingCache<SketchHitInfoSimple>,
    cache_out: &mut MappingCache<SketchHitInfoSimple>,
    index: &ReferenceIndex,
    poison_state: &mut PoisonState<'_>,
    strat: SkippingStrategy,
) -> bool
where
    Kmer<K>: KmerBits,
{
    poison_state.clear();
    cache_out.clear();

    poison_state.set_fragment_end(FragmentEnd::Left);
    let early_left =
        map_read::<K, SketchHitInfoSimple>(seq1, cache_left, hs, query, index, poison_state, strat);
    if poison_state.is_poisoned() {
        return false;
    }

    poison_state.set_fragment_end(FragmentEnd::Right);
    let early_right =
        map_read::<K, SketchHitInfoSimple>(seq2, cache_right, hs, query, index, poison_state, strat);
    if poison_state.is_poisoned() {
        return false;
    }

    merge_se_mappings(
        cache_left,
        cache_right,
        seq1.len() as i32,
        seq2.len() as i32,
        cache_out,
    );

    early_left || early_right
}

/// Map a single-end read fragment using bin-based hit collection (scATAC).
///
/// Returns `true` for early stop.
#[allow(clippy::too_many_arguments)]
pub fn map_se_fragment_atac<const K: usize>(
    seq: &[u8],
    hs: &mut HitSearcher<'_>,
    query: &mut PiscemStreamingQuery<'_, K>,
    cache_out: &mut MappingCache<SketchHitInfoSimple>,
    index: &ReferenceIndex,
    poison_state: &mut PoisonState<'_>,
    binning: &BinPos,
) -> bool
where
    Kmer<K>: KmerBits,
{
    poison_state.clear();
    map_read_atac::<K, SketchHitInfoSimple>(seq, cache_out, hs, query, index, poison_state, binning)
}

/// Map a paired-end read fragment using bin-based hit collection (scATAC).
///
/// Maps both ends independently with bin-based collection, then merges
/// using bin-aware pairing.
/// Returns `true` for early stop.
#[allow(clippy::too_many_arguments)]
pub fn map_pe_fragment_atac<const K: usize>(
    seq1: &[u8],
    seq2: &[u8],
    hs: &mut HitSearcher<'_>,
    query: &mut PiscemStreamingQuery<'_, K>,
    cache_left: &mut MappingCache<SketchHitInfoSimple>,
    cache_right: &mut MappingCache<SketchHitInfoSimple>,
    cache_out: &mut MappingCache<SketchHitInfoSimple>,
    index: &ReferenceIndex,
    poison_state: &mut PoisonState<'_>,
    binning: &BinPos,
) -> bool
where
    Kmer<K>: KmerBits,
{
    poison_state.clear();
    cache_out.clear();

    poison_state.set_fragment_end(FragmentEnd::Left);
    let early_left =
        map_read_atac::<K, SketchHitInfoSimple>(seq1, cache_left, hs, query, index, poison_state, binning);
    if poison_state.is_poisoned() {
        return false;
    }

    poison_state.set_fragment_end(FragmentEnd::Right);
    let early_right =
        map_read_atac::<K, SketchHitInfoSimple>(seq2, cache_right, hs, query, index, poison_state, binning);
    if poison_state.is_poisoned() {
        return false;
    }

    merge_se_mappings_binned(
        cache_left,
        cache_right,
        seq1.len() as i32,
        seq2.len() as i32,
        cache_out,
    );

    early_left || early_right
}
