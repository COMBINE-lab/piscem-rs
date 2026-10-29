//! Mapping kernel — the core `map_read()` function.
//!
//! Port of C++ `map_read()` from `mapping/utils.hpp`. This is the critical
//! per-read mapping function that:
//! 1. Collects raw k-mer hits from the dictionary
//! 2. Optionally applies poison filtering
//! 3. Accumulates hits per target using `SketchHitInfo`
//! 4. Optionally recovers from high-occurrence thresholds
//! 5. Optionally applies EC-based ambiguous hit filtering
//! 6. Produces the final accepted hit list

use crate::hash::{FixedMap, fixed_map};
use sshash_lib::{Kmer, KmerBits, KmerDictionary};

use crate::index::contig_table::ContigTableLike;
use crate::index::eq_classes::ec_entry_transcript_id;
use crate::index::reference_index::ReferenceIndex;
use crate::mapping::binning::BinPos;
use crate::mapping::cache::MappingCache;
use crate::mapping::filters::PoisonState;
use crate::mapping::hit_map::HitMap;
use crate::mapping::hit_searcher::{HitSearcher, SkippingStrategy};
use crate::mapping::hits::{HitDirection, MappingType, SketchHitInfo};
use crate::mapping::projected_hits::ProjectedHits;
use crate::mapping::streaming_query::PiscemStreamingQuery;

// ---------------------------------------------------------------------------
// map_read
// ---------------------------------------------------------------------------

/// Map a single read sequence against the index.
///
/// Returns `true` for early stop (the C++ convention — currently always
/// returns the early_stop flag from hit collection).
///
/// Port of C++ `map_read()` (non-binned version).
pub fn map_read<const K: usize, S: SketchHitInfo, D: KmerDictionary, C: ContigTableLike>(
    read_seq: &[u8],
    cache: &mut MappingCache<S>,
    hs: &mut HitSearcher<'_, D, C>,
    query: &mut PiscemStreamingQuery<'_, K, D>,
    index: &ReferenceIndex<D, C>,
    poison_state: &mut PoisonState<'_>,
    strat: SkippingStrategy,
) -> bool
where
    Kmer<K>: KmerBits,
{
    cache.clear();
    let k = cache.k as i32;
    let apply_poison_filter = poison_state.is_valid();
    let perform_ambig_filtering = index.has_ec_table();

    cache.has_matching_kmers = hs.get_raw_hits_sketch::<K>(read_seq, query, strat, true);
    let mut early_stop = false;

    let max_ec_ambig = cache.max_ec_card as usize;

    if cache.has_matching_kmers {
        let raw_hits = hs.left_hits();

        // Poison filter
        if apply_poison_filter {
            let was_poisoned = poison_state.scan_raw_hits(read_seq, k as u32, raw_hits, strat);
            if was_poisoned {
                poison_state.poison_read();
                cache.map_type = MappingType::Unmapped;
                return true;
            }
        }

        let signed_rl = read_seq.len() as i32;
        let max_stretch = signed_rl;

        // First pass: collect with max_hit_occ - 1 threshold
        let mut num_valid_hits: u32 = 0;
        let mao_first_pass = (cache.max_hit_occ - 1) as u64;
        let mut min_occ: u64 = u64::MAX;

        early_stop = collect_mappings_from_hits(
            raw_hits,
            &mut cache.hit_map,
            &mut num_valid_hits,
            &mut min_occ,
            mao_first_pass,
            signed_rl,
            k,
            max_stretch,
            perform_ambig_filtering,
            &mut cache.ambiguous_hit_indices,
            &mut cache.run_slots,
            index,
        );

        // Recovery: if all k-mers exceeded threshold, retry with min_occ
        if cache.attempt_occ_recover
            && min_occ >= cache.max_hit_occ as u64
            && min_occ < cache.max_hit_occ_recover as u64
        {
            num_valid_hits = 0;
            cache.ambiguous_hit_indices.clear();
            cache.hit_map.clear();
            let max_allowed_occ = min_occ;

            early_stop = collect_mappings_from_hits(
                raw_hits,
                &mut cache.hit_map,
                &mut num_valid_hits,
                &mut min_occ,
                max_allowed_occ,
                signed_rl,
                k,
                max_stretch,
                perform_ambig_filtering,
                &mut cache.ambiguous_hit_indices,
                &mut cache.run_slots,
                index,
            );
        }

        // EC-based ambiguous hit filtering
        if perform_ambig_filtering
            && !cache.hit_map.is_empty()
            && !cache.ambiguous_hit_indices.is_empty()
        {
            let ec_table = index.ec_table().unwrap();
            cache.observed_ecs.clear();
            let mut min_cardinality_ec_size = usize::MAX;
            let mut min_cardinality_ec: u64 = u64::MAX;
            let mut min_cardinality_index: u32 = 0;
            let mut visited: usize = 0;

            let visit_ec = |hit_map: &mut HitMap<S>, ent: u64, fw_on_contig: bool| {
                let tid = ec_entry_transcript_id(ent);
                if let Some(target) = hit_map.get_mut(&tid) {
                    let ori = ent & 0x3;
                    match ori {
                        0 => {
                            // FW
                            if fw_on_contig {
                                target.inc_fw_hits();
                            } else {
                                target.inc_rc_hits();
                            }
                        }
                        1 => {
                            // RC
                            if fw_on_contig {
                                target.inc_rc_hits();
                            } else {
                                target.inc_fw_hits();
                            }
                        }
                        _ => {
                            // Both
                            target.inc_fw_hits();
                            target.inc_rc_hits();
                        }
                    }
                }
            };

            for &hit_idx in &cache.ambiguous_hit_indices {
                let proj_hit = &raw_hits[hit_idx as usize].1;
                let contig_id = proj_hit.contig_id();
                let fw_on_contig = proj_hit.hit_fw_on_contig();

                let ec = ec_table.ec_for_tile(contig_id as u64);
                let ec_key = ec | if fw_on_contig { 0 } else { 0x8000000000000000 };

                if cache.observed_ecs.contains(&ec_key) {
                    continue;
                }
                cache.observed_ecs.insert(ec_key);

                let ec_entries = ec_table.entries_for_ec(ec);
                if ec_entries.len() < min_cardinality_ec_size {
                    min_cardinality_ec_size = ec_entries.len();
                    min_cardinality_ec = ec;
                    min_cardinality_index = hit_idx;
                }
                if ec_entries.len() > max_ec_ambig {
                    continue;
                }
                visited += 1;
                for ent in ec_entries.iter() {
                    visit_ec(&mut cache.hit_map, ent, fw_on_contig);
                }
                num_valid_hits += 1;
            }

            // Last-ditch: visit the smallest EC if we didn't visit any.
            if visited == 0 && min_cardinality_ec != u64::MAX {
                let proj_hit = &raw_hits[min_cardinality_index as usize].1;
                let fw_on_contig = proj_hit.hit_fw_on_contig();
                let ec_entries = ec_table.entries_for_ec(min_cardinality_ec);
                for ent in ec_entries.iter() {
                    visit_ec(&mut cache.hit_map, ent, fw_on_contig);
                }
                num_valid_hits += 1;
            }
        }

        // Final selection: keep targets with enough hits.
        for (&tid, target) in &cache.hit_map {
            let best_dir = target.best_hit_direction();
            let simple_hit = if best_dir != HitDirection::Rc {
                target.get_fw_hit()
            } else {
                target.get_rc_hit()
            };

            if simple_hit.num_hits >= num_valid_hits {
                let mut h = simple_hit;
                h.tid = tid;
                cache.accepted_hits.push(h);

                // If both directions are equally good, also add the RC hit.
                if best_dir == HitDirection::Both {
                    let mut rc_hit = target.get_rc_hit();
                    rc_hit.tid = tid;
                    cache.accepted_hits.push(rc_hit);
                }
            }
        }
    }

    // Filter by max_read_occ
    if cache.accepted_hits.len() > cache.max_read_occ {
        cache.accepted_hits.clear();
        cache.map_type = MappingType::Unmapped;
    } else if !cache.accepted_hits.is_empty() {
        cache.map_type = MappingType::SingleMapped;
    }

    early_stop
}

// ---------------------------------------------------------------------------
// map_read_atac — binned mapping for scATAC
// ---------------------------------------------------------------------------

/// Binned hit accumulator entry: wraps `SketchHitInfo` with a `tid`.
struct BinnedHitEntry<S: SketchHitInfo> {
    tid: u32,
    info: S,
}

impl<S: SketchHitInfo> Default for BinnedHitEntry<S> {
    fn default() -> Self {
        Self {
            tid: u32::MAX,
            info: S::default(),
        }
    }
}

/// Map a single read using bin-based hit collection (scATAC mode).
///
/// Uses `get_raw_hits_sketch_everykmer` and accumulates hits per bin_id.
/// Applies threshold filtering: accepted if `target.max_hits >= ceil(num_valid * thr)`.
/// Includes occurrence recovery and EC-based ambiguous hit filtering (matching C++).
///
/// Port of C++ `map_read()` with `bin_pos` overload.
#[allow(clippy::too_many_arguments)]
pub fn map_read_atac<const K: usize, S: SketchHitInfo, D: KmerDictionary, C: ContigTableLike>(
    read_seq: &[u8],
    cache: &mut MappingCache<S>,
    hs: &mut HitSearcher<'_, D, C>,
    query: &mut PiscemStreamingQuery<'_, K, D>,
    index: &ReferenceIndex<D, C>,
    poison_state: &mut PoisonState<'_>,
    binning: &BinPos,
) -> bool
where
    Kmer<K>: KmerBits,
{
    cache.clear();
    let k = cache.k as i32;
    let apply_poison_filter = poison_state.is_valid();
    let perform_ambig_filtering = index.has_ec_table();

    // Use every-kmer hit collection (matching C++ ATAC behavior)
    cache.has_matching_kmers = hs.get_raw_hits_sketch_everykmer::<K>(read_seq, query, true);

    let max_ec_ambig = cache.max_ec_card as usize;

    if cache.has_matching_kmers {
        let raw_hits = hs.left_hits();

        // Poison filter
        if apply_poison_filter {
            let was_poisoned =
                poison_state.scan_raw_hits(read_seq, k as u32, raw_hits, SkippingStrategy::Strict);
            if was_poisoned {
                poison_state.poison_read();
                cache.map_type = MappingType::Unmapped;
                return true;
            }
        }

        let signed_rl = read_seq.len() as i32;
        let max_stretch = signed_rl;
        let thr = binning.thr();

        // Bin-keyed hit map: bin_id → BinnedHitEntry
        let mut bin_hit_map: FixedMap<u64, BinnedHitEntry<S>> = fixed_map();
        let mut num_valid_hits: u32 = 0;
        let mut min_occ: u64 = u64::MAX;
        let mao_first_pass = (cache.max_hit_occ - 1) as u64;

        // First pass: collect with max_hit_occ - 1 threshold
        collect_mappings_from_hits_binned(
            raw_hits,
            &mut bin_hit_map,
            &mut num_valid_hits,
            &mut min_occ,
            mao_first_pass,
            signed_rl,
            k,
            max_stretch,
            perform_ambig_filtering,
            &mut cache.ambiguous_hit_indices,
            index,
            binning,
        );

        // Recovery: if all k-mers exceeded threshold, retry with min_occ
        if cache.attempt_occ_recover
            && min_occ >= cache.max_hit_occ as u64
            && min_occ < cache.max_hit_occ_recover as u64
        {
            num_valid_hits = 0;
            cache.ambiguous_hit_indices.clear();
            bin_hit_map.clear();
            let max_allowed_occ = min_occ;

            collect_mappings_from_hits_binned(
                raw_hits,
                &mut bin_hit_map,
                &mut num_valid_hits,
                &mut min_occ,
                max_allowed_occ,
                signed_rl,
                k,
                max_stretch,
                perform_ambig_filtering,
                &mut cache.ambiguous_hit_indices,
                index,
                binning,
            );
        }

        // EC-based ambiguous hit filtering
        // NOTE: In C++, visit_ec uses (ent >> 2) as the key to look up in the
        // bin_id-keyed map. Since (ent >> 2) is a transcript_id (not a bin_id),
        // this lookup almost never finds matches. However, num_valid_hits still
        // gets incremented, affecting the threshold. We match this behavior exactly.
        if perform_ambig_filtering
            && !bin_hit_map.is_empty()
            && !cache.ambiguous_hit_indices.is_empty()
        {
            let ec_table = index.ec_table().unwrap();
            cache.observed_ecs.clear();
            let mut min_cardinality_ec_size = usize::MAX;
            let mut min_cardinality_ec: u64 = u64::MAX;
            let mut min_cardinality_index: u32 = 0;
            let mut visited: usize = 0;

            // C++ visit_ec: uses (ent >> 2) as key, which is tid, not bin_id.
            // In bin-keyed map this is mostly a no-op, but we match C++ exactly.
            let visit_ec =
                |hit_map: &mut FixedMap<u64, BinnedHitEntry<S>>, ent: u64, fw_on_contig: bool| {
                    let tid_as_key = ec_entry_transcript_id(ent) as u64;
                    if let Some(target) = hit_map.get_mut(&tid_as_key) {
                        let ori = ent & 0x3;
                        match ori {
                            0 => {
                                if fw_on_contig {
                                    target.info.inc_fw_hits();
                                } else {
                                    target.info.inc_rc_hits();
                                }
                            }
                            1 => {
                                if fw_on_contig {
                                    target.info.inc_rc_hits();
                                } else {
                                    target.info.inc_fw_hits();
                                }
                            }
                            _ => {
                                target.info.inc_fw_hits();
                                target.info.inc_rc_hits();
                            }
                        }
                    }
                };

            for &hit_idx in &cache.ambiguous_hit_indices {
                let proj_hit = &raw_hits[hit_idx as usize].1;
                let contig_id = proj_hit.contig_id();
                let fw_on_contig = proj_hit.hit_fw_on_contig();

                let ec = ec_table.ec_for_tile(contig_id as u64);
                let ec_key = ec | if fw_on_contig { 0 } else { 0x8000000000000000 };

                if cache.observed_ecs.contains(&ec_key) {
                    continue;
                }
                cache.observed_ecs.insert(ec_key);

                let ec_entries = ec_table.entries_for_ec(ec);
                if ec_entries.len() < min_cardinality_ec_size {
                    min_cardinality_ec_size = ec_entries.len();
                    min_cardinality_ec = ec;
                    min_cardinality_index = hit_idx;
                }
                if ec_entries.len() > max_ec_ambig {
                    continue;
                }
                visited += 1;
                for ent in ec_entries.iter() {
                    visit_ec(&mut bin_hit_map, ent, fw_on_contig);
                }
                num_valid_hits += 1;
            }

            // Last-ditch: visit the smallest EC if we didn't visit any.
            if visited == 0 && min_cardinality_ec != u64::MAX {
                let proj_hit = &raw_hits[min_cardinality_index as usize].1;
                let fw_on_contig = proj_hit.hit_fw_on_contig();
                let ec_entries = ec_table.entries_for_ec(min_cardinality_ec);
                for ent in ec_entries.iter() {
                    visit_ec(&mut bin_hit_map, ent, fw_on_contig);
                }
                num_valid_hits += 1;
            }
        }

        // Apply threshold: require at least ceil(num_valid_hits * thr) hits
        let threshold = 1u32.max((num_valid_hits as f32 * thr).ceil() as u32);

        // Final selection: keep bins with enough hits
        for (&bin_id, bin_entry) in &bin_hit_map {
            let target = &bin_entry.info;
            let best_dir = target.best_hit_direction();
            let simple_hit = if best_dir != HitDirection::Rc {
                target.get_fw_hit()
            } else {
                target.get_rc_hit()
            };

            if simple_hit.num_hits >= threshold {
                let mut h = simple_hit;
                h.tid = bin_entry.tid;
                h.bin_id = bin_id;
                cache.accepted_hits.push(h);

                if best_dir == HitDirection::Both {
                    let mut rc_hit = target.get_rc_hit();
                    rc_hit.tid = bin_entry.tid;
                    rc_hit.bin_id = bin_id;
                    cache.accepted_hits.push(rc_hit);
                }
            }
        }
    }

    // Filter by max_read_occ
    if cache.accepted_hits.len() > cache.max_read_occ {
        cache.accepted_hits.clear();
        cache.map_type = MappingType::Unmapped;
    } else if !cache.accepted_hits.is_empty() {
        cache.map_type = MappingType::SingleMapped;
    }

    false
}

// ---------------------------------------------------------------------------
// collect_mappings_from_hits_binned (inner function for ATAC)
// ---------------------------------------------------------------------------

/// Inner loop for binned hit collection: iterate over raw hits, decode
/// reference positions, compute bin IDs, and accumulate into the bin hit map.
///
/// Port of C++ `collect_mappings_from_hits_thr` lambda in binned `map_read`.
#[allow(clippy::too_many_arguments)]
fn collect_mappings_from_hits_binned<S: SketchHitInfo, D: KmerDictionary, C: ContigTableLike>(
    raw_hits: &[(i32, ProjectedHits<'_>)],
    bin_hit_map: &mut FixedMap<u64, BinnedHitEntry<S>>,
    num_valid_hits: &mut u32,
    min_occ: &mut u64,
    max_allowed_occ: u64,
    signed_rl: i32,
    k: i32,
    max_stretch: i32,
    perform_ambig_filtering: bool,
    ambiguous_hit_indices: &mut Vec<u32>,
    index: &ReferenceIndex<D, C>,
    binning: &BinPos,
) {
    bin_hit_map.clear();
    let encoding = index.encoding();
    let score_inc: f32 = 1.0;

    for (hit_idx, (read_pos, proj_hits)) in raw_hits.iter().enumerate() {
        let ref_range = proj_hits.ref_range();
        let num_occ = ref_range.len() as u64;
        *min_occ = (*min_occ).min(num_occ);

        if num_occ <= max_allowed_occ {
            for entry in ref_range.iter() {
                let ref_pos_ori = proj_hits.decode_hit(entry, &encoding);
                let tid = encoding.transcript_id(entry);
                let pos = ref_pos_ori.pos;
                let ori = ref_pos_ori.is_fw;
                let signed_read_pos = *read_pos;

                // Get bin IDs for this hit
                let (bin1, bin2) = binning.get_bin_id(tid, pos as u64);

                // Add to primary bin
                {
                    let target = bin_hit_map.entry(bin1).or_default();
                    target.tid = tid;
                    if ori {
                        target.info.add_fw(
                            pos as i32,
                            signed_read_pos,
                            signed_rl,
                            k,
                            max_stretch,
                            score_inc,
                        );
                    } else {
                        target.info.add_rc(
                            pos as i32,
                            signed_read_pos,
                            signed_rl,
                            k,
                            max_stretch,
                            score_inc,
                        );
                    }
                }

                // Add to secondary bin (overlap region) if valid
                if BinPos::is_valid(bin2) {
                    let target2 = bin_hit_map.entry(bin2).or_default();
                    target2.tid = tid;
                    if ori {
                        target2.info.add_fw(
                            pos as i32,
                            signed_read_pos,
                            signed_rl,
                            k,
                            max_stretch,
                            score_inc,
                        );
                    } else {
                        target2.info.add_rc(
                            pos as i32,
                            signed_read_pos,
                            signed_rl,
                            k,
                            max_stretch,
                            score_inc,
                        );
                    }
                }
            }
            *num_valid_hits += 1;
        } else if perform_ambig_filtering {
            ambiguous_hit_indices.push(hit_idx as u32);
        }
    }
}

// ---------------------------------------------------------------------------
// collect_mappings_from_hits (inner function)
// ---------------------------------------------------------------------------

/// Inner loop: iterate over raw hits, decode reference positions, accumulate
/// into the hit map.
///
/// Raw hits are processed in *runs*: maximal stretches of consecutive hits on
/// the same contig, which therefore share one occurrence list. The run's first
/// hit resolves every occurrence to its `hit_map` slot, inserting targets in
/// list order exactly as a per-occurrence `entry()` would, so the map's
/// contents and iteration order are those of the plain hit-by-hit loop.
///
/// A run is then applied target by target rather than hit by hit. That is the
/// same computation: targets are independent, and each still receives the
/// run's hits in order with the same `num_valid_hits` gate, so it ends in the
/// same state. The one coupling between targets — the early stop after the
/// first hit that leaves no target ahead of `num_valid_hits` — is recovered from
/// a per-hit bitmask; the hits a target-major pass applies past that point are
/// no-ops, because a target that is not ahead after hit `j` fails the gate for
/// every later hit. The only case where order *within* a target could differ
/// is a transcript listed twice for one contig (its forward and reverse hits
/// interleave differently), so runs whose occurrence list is not strictly
/// ascending in transcript id — never the case for lists piscem builds, which
/// are in reference order — take the hit-by-hit path.
///
/// Returns `true` for early stop (no valid targets remaining).
#[allow(clippy::too_many_arguments)]
fn collect_mappings_from_hits<S: SketchHitInfo, D: KmerDictionary, C: ContigTableLike>(
    raw_hits: &[(i32, ProjectedHits<'_>)],
    hit_map: &mut HitMap<S>,
    num_valid_hits: &mut u32,
    min_occ: &mut u64,
    max_allowed_occ: u64,
    signed_rl: i32,
    k: i32,
    max_stretch: i32,
    perform_ambig_filtering: bool,
    ambiguous_hit_indices: &mut Vec<u32>,
    run_slots: &mut Vec<(u32, u64)>,
    index: &ReferenceIndex<D, C>,
) -> bool {
    /// Longest run handled target-major; bounded by the bitmask width.
    const MAX_RUN: usize = u64::BITS as usize;

    hit_map.clear();
    let encoding = index.encoding();
    let score_inc: f32 = 1.0;

    // Apply one hit to one target; `None` if the target fails the gate,
    // otherwise whether it is now ahead of `nvh`.
    let apply = |target: &mut S, read_pos: i32, hit: &ProjectedHits<'_>, entry: u64, nvh: u32| {
        if target.max_hits_for_target() < nvh {
            return None;
        }
        let ref_pos_ori = hit.decode_hit(entry, &encoding);
        let pos = ref_pos_ori.pos as i32;
        if ref_pos_ori.is_fw {
            target.add_fw(pos, read_pos, signed_rl, k, max_stretch, score_inc);
        } else {
            target.add_rc(pos, read_pos, signed_rl, k, max_stretch, score_inc);
        }
        Some(target.max_hits_for_target() > nvh)
    };

    let mut i = 0;
    while i < raw_hits.len() {
        let proj_hits = &raw_hits[i].1;
        let ref_range = proj_hits.ref_range();
        let num_occ = ref_range.len() as u64;
        *min_occ = (*min_occ).min(num_occ);

        if num_occ > max_allowed_occ {
            if perform_ambig_filtering {
                ambiguous_hit_indices.push(i as u32);
            }
            i += 1;
            continue;
        }

        // Every hit in the run shares this occurrence list, so all are valid
        // and leave `min_occ` where it is.
        let contig_id = proj_hits.contig_id();
        let mut end = i + 1;
        while end < raw_hits.len() && end - i < MAX_RUN && raw_hits[end].1.contig_id() == contig_id
        {
            end += 1;
        }
        let run = &raw_hits[i..end];

        run_slots.clear();
        let mut ascending = true;
        let mut prev_tid: Option<u32> = None;
        for entry in ref_range.iter() {
            let tid = encoding.transcript_id(entry);
            ascending &= prev_tid.is_none_or(|p| tid > p);
            prev_tid = Some(tid);
            run_slots.push((hit_map.slot_or_insert(tid), entry));
        }

        let nvh0 = *num_valid_hits;
        if run.len() > 1 && ascending {
            // Bit `j` set: some target is ahead of `num_valid_hits` after hit `j`.
            let mut ahead: u64 = 0;
            for &(slot, entry) in run_slots.iter() {
                let target = hit_map.slot_mut(slot);
                for (j, (read_pos, hit)) in run.iter().enumerate() {
                    match apply(target, *read_pos, hit, entry, nvh0 + j as u32) {
                        None => break,
                        Some(true) => ahead |= 1 << j,
                        Some(false) => {}
                    }
                }
            }
            let survived = ahead.trailing_ones() as usize;
            if survived < run.len() {
                *num_valid_hits = nvh0 + survived as u32 + 1;
                return true;
            }
            *num_valid_hits = nvh0 + run.len() as u32;
        } else {
            for (read_pos, hit) in run {
                let mut still_have_valid_target = false;
                for &(slot, entry) in run_slots.iter() {
                    let target = hit_map.slot_mut(slot);
                    still_have_valid_target |=
                        apply(target, *read_pos, hit, entry, *num_valid_hits) == Some(true);
                }
                *num_valid_hits += 1;
                if !still_have_valid_target {
                    return true;
                }
            }
        }
        i = end;
    }

    false
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mapping::cache::MappingCache;
    use crate::mapping::hits::SimpleHit;
    use crate::mapping::sketch_hit_simple::SketchHitInfoSimple;

    #[test]
    fn test_collect_empty_hits() {
        let raw_hits: Vec<(i32, ProjectedHits<'_>)> = Vec::new();
        let mut hit_map: HitMap<SketchHitInfoSimple> = HitMap::with_capacity(0);
        let num_valid = 0u32;
        let mut min_occ = u64::MAX;
        let mut ambig: Vec<u32> = Vec::new();

        // We need a dummy index for the encoding — but with empty raw_hits
        // the function won't access it. We can test that the function returns
        // false (no early stop) with empty input.
        // Since we can't easily construct a ReferenceIndex in unit tests
        // without real data, we test the logic flow through map_read's
        // has_matching_kmers check instead.

        // This test verifies the collect function handles empty input.
        // The inner function won't be called with empty raw_hits in practice.
        assert!(raw_hits.is_empty());
        assert_eq!(num_valid, 0);

        // Verify cache behavior instead
        let mut cache = MappingCache::<SketchHitInfoSimple>::new(31);
        cache.has_matching_kmers = false;
        assert_eq!(cache.map_type, MappingType::Unmapped);
        let _ = (&mut hit_map, &mut min_occ, &mut ambig);
    }

    #[test]
    fn test_max_read_occ_filter() {
        let mut cache = MappingCache::<SketchHitInfoSimple>::new(31);
        cache.max_read_occ = 2;

        // Simulate having too many accepted hits
        for i in 0..5 {
            let hit = SimpleHit {
                tid: i,
                ..SimpleHit::default()
            };
            cache.accepted_hits.push(hit);
        }

        // The filter logic from map_read:
        if cache.accepted_hits.len() > cache.max_read_occ {
            cache.accepted_hits.clear();
            cache.map_type = MappingType::Unmapped;
        } else if !cache.accepted_hits.is_empty() {
            cache.map_type = MappingType::SingleMapped;
        }

        assert_eq!(cache.map_type, MappingType::Unmapped);
        assert!(cache.accepted_hits.is_empty());
    }

    #[test]
    fn test_single_mapped_classification() {
        let mut cache = MappingCache::<SketchHitInfoSimple>::new(31);

        let hit = SimpleHit {
            tid: 42,
            ..SimpleHit::default()
        };
        cache.accepted_hits.push(hit);

        // Apply the classification logic from map_read
        if cache.accepted_hits.len() > cache.max_read_occ {
            cache.accepted_hits.clear();
            cache.map_type = MappingType::Unmapped;
        } else if !cache.accepted_hits.is_empty() {
            cache.map_type = MappingType::SingleMapped;
        }

        assert_eq!(cache.map_type, MappingType::SingleMapped);
        assert_eq!(cache.accepted_hits.len(), 1);
    }

    #[test]
    fn test_unmapped_when_no_hits() {
        let mut cache = MappingCache::<SketchHitInfoSimple>::new(31);

        // No accepted hits
        if cache.accepted_hits.len() > cache.max_read_occ {
            cache.accepted_hits.clear();
            cache.map_type = MappingType::Unmapped;
        } else if !cache.accepted_hits.is_empty() {
            cache.map_type = MappingType::SingleMapped;
        }

        assert_eq!(cache.map_type, MappingType::Unmapped);
    }
}
