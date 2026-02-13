//! Paired-end merge — combines left and right single-end mappings.
//!
//! Port of C++ `merge_se_mappings()` (non-binned version). Looks for
//! concordant read pairs (left-FW × right-RC and left-RC × right-FW)
//! on the same transcript, with fragment length constraints.

use crate::mapping::cache::MappingCache;
use crate::mapping::hits::{MappingType, SimpleHit, SketchHitInfo};

// ---------------------------------------------------------------------------
// merge_se_mappings
// ---------------------------------------------------------------------------

/// Merge left and right single-end mapping results into paired-end output.
///
/// Algorithm:
/// 1. Sort both `accepted_hits` by `(orientation desc, tid asc, pos asc)`
/// 2. Split into FW and RC sublists via partition point
/// 3. Two-pointer merge: left_FW × right_RC, left_RC × right_FW
/// 4. Fragment length check: `-32 < frag_len < 2000`
/// 5. Output: `MappedPair`, `MappedFirstOrphan`, `MappedSecondOrphan`, or `Unmapped`
///
/// Port of C++ `merge_se_mappings()`.
pub fn merge_se_mappings<S: SketchHitInfo>(
    cache_left: &mut MappingCache<S>,
    cache_right: &mut MappingCache<S>,
    left_len: i32,
    right_len: i32,
    cache_out: &mut MappingCache<S>,
) {
    cache_out.clear();

    let had_matching_kmers_left = cache_left.has_matching_kmers;
    let had_matching_kmers_right = cache_right.has_matching_kmers;

    let num_accepted_left = cache_left.accepted_hits.len();
    let num_accepted_right = cache_right.accepted_hits.len();

    if num_accepted_left > 0 && num_accepted_right > 0 {
        // Sort by (orientation desc, tid asc, pos asc)
        cache_left.accepted_hits.sort_by(simple_hit_cmp);
        cache_right.accepted_hits.sort_by(simple_hit_cmp);

        // Sentinel = smallest possible RC element in sort order.
        // Sort: FW first (is_fw desc), then tid asc, then pos asc.
        // All FW elements compare Less than this; all RC elements compare >= this.
        let sentinel = SimpleHit {
            is_fw: false,
            tid: 0,
            pos: i32::MIN,
            ..SimpleHit::default()
        };

        // Split left into FW and RC sublists
        let left_fw_end = cache_left
            .accepted_hits
            .partition_point(|h| simple_hit_cmp(h, &sentinel) == std::cmp::Ordering::Less);
        let (left_fw, left_rc) = cache_left.accepted_hits.split_at(left_fw_end);

        // Split right into FW and RC sublists
        let right_fw_end = cache_right
            .accepted_hits
            .partition_point(|h| simple_hit_cmp(h, &sentinel) == std::cmp::Ordering::Less);
        let (right_fw, right_rc) = cache_right.accepted_hits.split_at(right_fw_end);



        // Merge left_FW × right_RC
        merge_lists(
            left_fw,
            right_rc,
            left_len,
            right_len,
            &mut cache_out.accepted_hits,
        );

        // Merge left_RC × right_FW
        merge_lists(
            left_rc,
            right_fw,
            left_len,
            right_len,
            &mut cache_out.accepted_hits,
        );

        cache_out.map_type = if !cache_out.accepted_hits.is_empty() {
            MappingType::MappedPair
        } else {
            MappingType::Unmapped
        };
    } else if num_accepted_left > 0 && !had_matching_kmers_right {
        std::mem::swap(
            &mut cache_left.accepted_hits,
            &mut cache_out.accepted_hits,
        );
        cache_out.map_type = if !cache_out.accepted_hits.is_empty() {
            MappingType::MappedFirstOrphan
        } else {
            MappingType::Unmapped
        };
    } else if num_accepted_right > 0 && !had_matching_kmers_left {
        std::mem::swap(
            &mut cache_right.accepted_hits,
            &mut cache_out.accepted_hits,
        );
        cache_out.map_type = if !cache_out.accepted_hits.is_empty() {
            MappingType::MappedSecondOrphan
        } else {
            MappingType::Unmapped
        };
    }
    // else: both empty or both had matching k-mers but no accepted → Unmapped
}

/// Sort comparator for `SimpleHit`: orientation desc, tid asc, pos asc.
fn simple_hit_cmp(a: &SimpleHit, b: &SimpleHit) -> std::cmp::Ordering {
    // FW sorts before RC (is_fw true > false, so reverse)
    b.is_fw
        .cmp(&a.is_fw)
        .then(a.tid.cmp(&b.tid))
        .then(a.pos.cmp(&b.pos))
}

/// Two-pointer merge of compatible orientation lists.
///
/// `list1` hits have one orientation, `list2` hits have the opposite.
/// Looks for pairs on the same tid with acceptable fragment length.
fn merge_lists(
    list1: &[SimpleHit],
    list2: &[SimpleHit],
    left_len: i32,
    right_len: i32,
    out: &mut Vec<SimpleHit>,
) {
    let mut i1 = 0;
    let mut i2 = 0;

    while i1 < list1.len() && i2 < list2.len() {
        if list1[i1].tid < list2[i2].tid {
            i1 += 1;
        } else {
            if list1[i1].tid == list2[i2].tid {
                let pos_fw = if list1[i1].is_fw {
                    list1[i1].pos
                } else {
                    list2[i2].pos
                };
                let pos_rc = if list1[i1].is_fw {
                    list2[i2].pos
                } else {
                    list1[i1].pos
                };
                let frag_len = pos_rc - pos_fw;

                if frag_len > -32 && frag_len < 2000 {
                    let right_is_rc = !list2[i2].is_fw;
                    let tlen = if right_is_rc {
                        (list2[i2].pos + right_len - list1[i1].pos) + 1
                    } else {
                        (list1[i1].pos + left_len - list2[i2].pos) + 1
                    };
                    out.push(SimpleHit {
                        is_fw: list1[i1].is_fw,
                        mate_is_fw: list2[i2].is_fw,
                        pos: list1[i1].pos,
                        score: 0.0,
                        num_hits: 0,
                        tid: list1[i1].tid,
                        mate_pos: list2[i2].pos,
                        fragment_length: tlen,
                    });
                    i1 += 1;
                }
            }
            i2 += 1;
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mapping::sketch_hit_simple::SketchHitInfoSimple;

    fn make_hit(tid: u32, pos: i32, is_fw: bool) -> SimpleHit {
        SimpleHit {
            is_fw,
            mate_is_fw: false,
            pos,
            score: 1.0,
            num_hits: 1,
            tid,
            ..SimpleHit::default()
        }
    }

    #[test]
    fn test_merge_both_empty() {
        let mut left = MappingCache::<SketchHitInfoSimple>::new(31);
        let mut right = MappingCache::<SketchHitInfoSimple>::new(31);
        let mut out = MappingCache::<SketchHitInfoSimple>::new(31);

        merge_se_mappings(&mut left, &mut right, 100, 100, &mut out);
        assert_eq!(out.map_type, MappingType::Unmapped);
        assert!(out.accepted_hits.is_empty());
    }

    #[test]
    fn test_merge_left_orphan() {
        let mut left = MappingCache::<SketchHitInfoSimple>::new(31);
        let mut right = MappingCache::<SketchHitInfoSimple>::new(31);
        let mut out = MappingCache::<SketchHitInfoSimple>::new(31);

        left.accepted_hits.push(make_hit(0, 100, true));
        // Right had no matching k-mers at all
        right.has_matching_kmers = false;

        merge_se_mappings(&mut left, &mut right, 100, 100, &mut out);
        assert_eq!(out.map_type, MappingType::MappedFirstOrphan);
        assert_eq!(out.accepted_hits.len(), 1);
    }

    #[test]
    fn test_merge_right_orphan() {
        let mut left = MappingCache::<SketchHitInfoSimple>::new(31);
        let mut right = MappingCache::<SketchHitInfoSimple>::new(31);
        let mut out = MappingCache::<SketchHitInfoSimple>::new(31);

        left.has_matching_kmers = false;
        right.accepted_hits.push(make_hit(0, 200, false));

        merge_se_mappings(&mut left, &mut right, 100, 100, &mut out);
        assert_eq!(out.map_type, MappingType::MappedSecondOrphan);
        assert_eq!(out.accepted_hits.len(), 1);
    }

    #[test]
    fn test_merge_concordant_pair() {
        let mut left = MappingCache::<SketchHitInfoSimple>::new(31);
        let mut right = MappingCache::<SketchHitInfoSimple>::new(31);
        let mut out = MappingCache::<SketchHitInfoSimple>::new(31);

        // Left FW at pos 100, Right RC at pos 300 on same tid
        left.accepted_hits.push(make_hit(5, 100, true));
        left.has_matching_kmers = true;
        right.accepted_hits.push(make_hit(5, 300, false));
        right.has_matching_kmers = true;

        merge_se_mappings(&mut left, &mut right, 100, 100, &mut out);
        assert_eq!(out.map_type, MappingType::MappedPair);
        assert_eq!(out.accepted_hits.len(), 1);
        assert_eq!(out.accepted_hits[0].tid, 5);
        assert!(out.accepted_hits[0].is_fw);
        assert!(!out.accepted_hits[0].mate_is_fw);
    }

    #[test]
    fn test_merge_frag_len_filter() {
        let mut left = MappingCache::<SketchHitInfoSimple>::new(31);
        let mut right = MappingCache::<SketchHitInfoSimple>::new(31);
        let mut out = MappingCache::<SketchHitInfoSimple>::new(31);

        // Fragment length too large (> 2000)
        left.accepted_hits.push(make_hit(5, 100, true));
        left.has_matching_kmers = true;
        right.accepted_hits.push(make_hit(5, 3000, false));
        right.has_matching_kmers = true;

        merge_se_mappings(&mut left, &mut right, 100, 100, &mut out);
        // frag_len = 3000 - 100 = 2900 > 2000 → rejected
        assert_eq!(out.map_type, MappingType::Unmapped);
        assert!(out.accepted_hits.is_empty());
    }

    #[test]
    fn test_merge_different_tids() {
        let mut left = MappingCache::<SketchHitInfoSimple>::new(31);
        let mut right = MappingCache::<SketchHitInfoSimple>::new(31);
        let mut out = MappingCache::<SketchHitInfoSimple>::new(31);

        left.accepted_hits.push(make_hit(1, 100, true));
        left.has_matching_kmers = true;
        right.accepted_hits.push(make_hit(2, 300, false));
        right.has_matching_kmers = true;

        merge_se_mappings(&mut left, &mut right, 100, 100, &mut out);
        // Different tids → no pairs
        assert_eq!(out.map_type, MappingType::Unmapped);
    }
}
