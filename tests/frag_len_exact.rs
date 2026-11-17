//! Fragment length and position conventions, end to end on a real index.
//!
//! Builds a one-transcript index from a hand-written cDBG (the transcript is a
//! single unitig), maps read pairs cut from known coordinates, and pins the
//! exact `pos`, `mate_pos` and `fragment_length` the mapper reports, plus the
//! start and length the bulk and scATAC RAD writers emit.
//!
//! Convention: positions are 0-based leftmost reference coordinates, and a
//! fragment covering the half-open interval `[start, end)` has length
//! `end - start`.

use std::io::Write;
use std::path::Path;

use piscem_rs::index::build::{BuildConfig, build_index};
use piscem_rs::index::reference_index::ReferenceIndex;
use piscem_rs::io::rad::{RadWriter, write_atac_record, write_bulk_record};
use piscem_rs::mapping::binning::BinPos;
use piscem_rs::mapping::cache::MappingCache;
use piscem_rs::mapping::filters::PoisonState;
use piscem_rs::mapping::hit_searcher::{HitSearcher, SkippingStrategy};
use piscem_rs::mapping::hits::{MappingType, SimpleHit};
use piscem_rs::mapping::map_fragment::{
    map_pe_fragment, map_pe_fragment_atac, map_se_fragment_atac,
};
use piscem_rs::mapping::merge_pairs::{remove_duplicate_hits_pub, simple_hit_cmp_bins};
use piscem_rs::mapping::overlap::{OverlapType, find_overlap, merged_mate_pos};
use piscem_rs::mapping::sketch_hit_simple::SketchHitInfoSimple;
use piscem_rs::mapping::streaming_query::PiscemStreamingQuery;

const K: usize = 31;
const TXP_LEN: usize = 2000;

/// Deterministic high-complexity DNA (LCG over the 4 bases).
fn gen_seq(n: usize, seed: u64) -> Vec<u8> {
    const B: [u8; 4] = *b"ACGT";
    let mut x = seed;
    (0..n)
        .map(|_| {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            B[((x >> 33) & 3) as usize]
        })
        .collect()
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            o => o,
        })
        .collect()
}

/// Write a cDBG for `seq` as one unitig and build a piscem index from it.
fn build_one_txp_index(dir: &Path, seq: &[u8]) -> ReferenceIndex {
    let cf = dir.join("cf");
    let mut f = std::fs::File::create(cf.with_extension("cf_seg")).unwrap();
    write!(f, "1\t").unwrap();
    f.write_all(seq).unwrap();
    writeln!(f).unwrap();
    std::fs::write(
        cf.with_extension("cf_seq"),
        "Reference:1_Sequence:tx0\t1+\n",
    )
    .unwrap();
    std::fs::write(cf.with_extension("json"), "{}").unwrap();
    let out = dir.join("idx").join("tx");
    std::fs::create_dir_all(out.parent().unwrap()).unwrap();
    let config = BuildConfig {
        input_prefix: cf,
        output_prefix: out.clone(),
        k: K,
        m: 19,
        build_ec_table: false,
        num_threads: 1,
        seed: 1,
        mphf_build: Default::default(),
        emit_tiny: Some(false),
        tmp_dir: Some(dir.join("sshash_tmp")),
        ram_limit_gib: Some(1),
    };
    build_index(&config).expect("build index");
    ReferenceIndex::load(&out, false, false).expect("load index")
}

struct Fixture {
    _dir: tempfile::TempDir,
    txp: Vec<u8>,
    index: ReferenceIndex,
}

fn fixture() -> Fixture {
    let dir = tempfile::tempdir().unwrap();
    let txp = gen_seq(TXP_LEN, 0xC0FFEE);
    let index = build_one_txp_index(dir.path(), &txp);
    assert_eq!(index.num_refs(), 1);
    assert_eq!(index.ref_len(0), TXP_LEN as u64);
    Fixture {
        _dir: dir,
        txp,
        index,
    }
}

/// Map a pair through the bulk (`merge_se_mappings`) path.
fn map_bulk(index: &ReferenceIndex, r1: &[u8], r2: &[u8]) -> (MappingType, Vec<SimpleHit>) {
    let mut hs = HitSearcher::new(index);
    let mut query = PiscemStreamingQuery::<K, _>::new(index.dict());
    let mut left = MappingCache::<SketchHitInfoSimple>::new(K);
    let mut right = MappingCache::<SketchHitInfoSimple>::new(K);
    let mut out = MappingCache::<SketchHitInfoSimple>::new(K);
    let mut poison = PoisonState::new(index.poison_table());
    poison.paired_for_mapping = true;
    map_pe_fragment(
        r1,
        r2,
        &mut hs,
        &mut query,
        &mut left,
        &mut right,
        &mut out,
        index,
        &mut poison,
        SkippingStrategy::Permissive,
    );
    (out.map_type, out.accepted_hits.clone())
}

/// Map a pair through the scATAC (`merge_se_mappings_binned`) path.
fn map_atac(index: &ReferenceIndex, r1: &[u8], r2: &[u8]) -> (MappingType, Vec<SimpleHit>) {
    let binning = BinPos::new(index, 1000, 300, 0.7);
    let mut hs = HitSearcher::new(index);
    let mut query = PiscemStreamingQuery::<K, _>::new(index.dict());
    let mut left = MappingCache::<SketchHitInfoSimple>::new(K);
    let mut right = MappingCache::<SketchHitInfoSimple>::new(K);
    let mut out = MappingCache::<SketchHitInfoSimple>::new(K);
    let mut poison = PoisonState::new(index.poison_table());
    poison.paired_for_mapping = true;
    map_pe_fragment_atac(
        r1,
        r2,
        &mut hs,
        &mut query,
        &mut left,
        &mut right,
        &mut out,
        index,
        &mut poison,
        &binning,
    );
    (out.map_type, out.accepted_hits.clone())
}

/// Map an overlapping pair the way `map-scatac` does: merge the mates, map
/// the merged fragment as one read, then fill in the pair fields (this mirrors
/// the overlap branch of the scATAC processor).
fn map_atac_overlap(index: &ReferenceIndex, r1: &[u8], r2: &[u8]) -> (OverlapType, Vec<SimpleHit>) {
    let binning = BinPos::new(index, 1000, 300, 0.7);
    let mut hs = HitSearcher::new(index);
    let mut query = PiscemStreamingQuery::<K, _>::new(index.dict());
    let mut out = MappingCache::<SketchHitInfoSimple>::new(K);
    let mut poison = PoisonState::new(index.poison_table());
    let mate_ov = find_overlap(r1, r2, 30, 0);
    assert_ne!(mate_ov.ov_type, OverlapType::NoOverlap);
    poison.paired_for_mapping = false;
    map_se_fragment_atac(
        &mate_ov.frag,
        &mut hs,
        &mut query,
        &mut out,
        index,
        &mut poison,
        &binning,
    );
    out.accepted_hits.sort_by(simple_hit_cmp_bins);
    remove_duplicate_hits_pub(&mut out.accepted_hits);
    let longer_len = r1.len().max(r2.len()) as i32;
    for hit in &mut out.accepted_hits {
        hit.fragment_length = mate_ov.frag_length as i32;
        hit.mate_is_fw = !hit.is_fw;
        hit.mate_pos = merged_mate_pos(hit.pos, hit.is_fw, hit.fragment_length, longer_len);
        if mate_ov.ov_type == OverlapType::Dovetail {
            hit.mate_pos = hit.pos;
        }
    }
    (mate_ov.ov_type, out.accepted_hits)
}

/// (start, frag_len) of each alignment of a bulk RAD pair record.
fn bulk_rad(hits: &[SimpleHit]) -> Vec<(u32, u16)> {
    let mut w = RadWriter::new();
    write_bulk_record(MappingType::MappedPair, hits, &mut w);
    w.as_bytes()[5..]
        .as_chunks::<10>()
        .0
        .iter()
        .map(|a| {
            (
                u32::from_le_bytes(a[4..8].try_into().unwrap()),
                u16::from_le_bytes(a[8..10].try_into().unwrap()),
            )
        })
        .collect()
}

/// (start, frag_len) of each alignment of a scATAC RAD pair record.
fn atac_rad(hits: &[SimpleHit], tn5_shift: bool) -> Vec<(u32, u16)> {
    let mut w = RadWriter::new();
    write_atac_record(0, 16, MappingType::MappedPair, hits, tn5_shift, &mut w);
    w.as_bytes()[8..]
        .as_chunks::<11>()
        .0
        .iter()
        .map(|a| {
            (
                u32::from_le_bytes(a[5..9].try_into().unwrap()),
                u16::from_le_bytes(a[9..11].try_into().unwrap()),
            )
        })
        .collect()
}

/// One pair cut from `[50, 430)`: r1 = t[50..130] forward and r2 = the reverse
/// complement of t[350..430]. Hit positions are exact 0-based leftmost
/// coordinates, the fragment is 380 bases (not 381), and both read orders and
/// both merge paths agree.
#[test]
fn pair_fragment_length_is_exact() {
    let fx = fixture();
    let t = &fx.txp;
    let fw = &t[50..130];
    let rc = revcomp(&t[350..430]);
    for (name, map) in [
        ("bulk", map_bulk as fn(&ReferenceIndex, &[u8], &[u8]) -> _),
        ("atac", map_atac),
    ] {
        // r1 forward, r2 reverse complement
        let (mt, hits) = map(&fx.index, fw, &rc);
        assert_eq!(mt, MappingType::MappedPair, "{name}");
        assert_eq!(hits.len(), 1, "{name}");
        let h = &hits[0];
        assert_eq!(
            (h.is_fw, h.pos, h.mate_is_fw, h.mate_pos),
            (true, 50, false, 350),
            "{name}"
        );
        assert_eq!(h.fragment_length, 380, "{name} {h:?}");
        // r1 reverse complement, r2 forward
        let (mt, hits) = map(&fx.index, &rc, fw);
        assert_eq!(mt, MappingType::MappedPair, "{name}");
        let h = &hits[0];
        assert_eq!(
            (h.is_fw, h.pos, h.mate_is_fw, h.mate_pos),
            (false, 350, true, 50),
            "{name}"
        );
        assert_eq!(h.fragment_length, 380, "{name} {h:?}");
    }

    // RAD: bulk writes start 50, length 380; scATAC with the Tn5 shift writes
    // start 54 and length 371, i.e. the fragments.tsv interval [54, 425)
    // (+4 on the forward read's start, -5 on the reverse read's end).
    let (_, hits) = map_bulk(&fx.index, fw, &rc);
    assert_eq!(bulk_rad(&hits), [(50, 380)]);
    let (_, hits) = map_atac(&fx.index, &rc, fw);
    assert_eq!(atac_rad(&hits, false), [(50, 380)]);
    assert_eq!(atac_rad(&hits, true), [(54, 371)]);
}

/// A forward mate that starts 10 bases before the transcript (`pos = -10`):
/// the pair spans [-10, 380), 390 bases, and the RAD writer clips it to
/// [0, 380), 380 bases.
#[test]
fn overhanging_pair_is_clipped_exactly() {
    let fx = fixture();
    let t = &fx.txp;
    let mut r1 = gen_seq(10, 7);
    r1.extend_from_slice(&t[0..70]);
    let r2 = revcomp(&t[300..380]);
    let (mt, hits) = map_bulk(&fx.index, &r1, &r2);
    assert_eq!(mt, MappingType::MappedPair);
    assert_eq!((hits[0].pos, hits[0].mate_pos), (-10, 300));
    assert_eq!(hits[0].fragment_length, 390);
    assert_eq!(bulk_rad(&hits), [(0, 380)]);
}

/// scATAC mates that overlap are merged into one fragment and mapped as one
/// read. The fragment [200, 380) must be written as start 200, length 180,
/// whichever strand it maps to (the reverse-complement case used to start at
/// 139, `frag_len - shorter_len + 1` bases too far left).
#[test]
fn overlapping_atac_pair_is_exact() {
    let fx = fixture();
    let t = &fx.txp;
    let a = t[200..320].to_vec();
    let b = revcomp(&t[260..380]);
    for (name, r1, r2) in [("fw", &a, &b), ("rc", &b, &a)] {
        let (ov, hits) = map_atac_overlap(&fx.index, r1, r2);
        assert_eq!(ov, OverlapType::Overlap, "{name}");
        assert_eq!(hits.len(), 1, "{name}");
        assert_eq!(hits[0].is_fw, name == "fw");
        assert_eq!(hits[0].pos, 200, "{name}");
        assert_eq!(hits[0].fragment_length, 180, "{name}");
        assert_eq!(atac_rad(&hits, false), [(200, 180)], "{name} {hits:?}");
        assert_eq!(atac_rad(&hits, true), [(204, 171)], "{name}");
    }
    // Unequal mates: a 100-base r1 and a 130-base r2 over [200, 400).
    let a = t[200..300].to_vec();
    let b = revcomp(&t[270..400]);
    for (name, r1, r2) in [("fw", &a, &b), ("rc", &b, &a)] {
        let (ov, hits) = map_atac_overlap(&fx.index, r1, r2);
        assert_eq!(ov, OverlapType::Overlap, "{name}");
        assert_eq!(atac_rad(&hits, false), [(200, 200)], "{name} {hits:?}");
    }
    // Read-through (dovetail): an 80-base fragment read past its ends by both
    // 100-base mates is trimmed to [200, 280).
    let mut a = t[200..280].to_vec();
    a.extend(gen_seq(20, 3));
    let mut b = revcomp(&t[200..280]);
    b.extend(gen_seq(20, 4));
    let (ov, hits) = map_atac_overlap(&fx.index, &a, &b);
    assert_eq!(ov, OverlapType::Dovetail);
    assert_eq!(atac_rad(&hits, false), [(200, 80)]);
}

/// Adapter-like bases that match nowhere in the transcript.
fn adapter(n: usize, seed: u64) -> Vec<u8> {
    gen_seq(n, seed ^ 0xADA9_7E55)
}

/// Dovetailed pairs: the 80-base insert [200, 280) read by 95-base mates that
/// run 15 bases into adapter (the scATAC merge pairs mates whose starts differ
/// by less than 20). The reverse mate's projected start (185) lies before the
/// forward mate's (200); the fragment is still [200, 280), the span between
/// the two 5' ends, in both read orders and on both merge paths.
#[test]
fn dovetailed_pair_starts_at_forward_mate() {
    let fx = fixture();
    let t = &fx.txp;
    let mut fw = t[200..280].to_vec();
    fw.extend(adapter(15, 1));
    let mut rc = revcomp(&t[200..280]);
    rc.extend(adapter(15, 2));
    for (name, map) in [
        ("bulk", map_bulk as fn(&ReferenceIndex, &[u8], &[u8]) -> _),
        ("atac", map_atac),
    ] {
        for (order, r1, r2) in [("fw/rc", &fw, &rc), ("rc/fw", &rc, &fw)] {
            let (mt, hits) = map(&fx.index, r1, r2);
            assert_eq!(mt, MappingType::MappedPair, "{name} {order}");
            assert_eq!(hits.len(), 1, "{name} {order}");
            let h = &hits[0];
            assert_eq!(h.pair_start(), 200, "{name} {order} {h:?}");
            assert_eq!(h.fragment_length, 80, "{name} {order} {h:?}");
            if name == "bulk" {
                assert_eq!(bulk_rad(&hits), [(200, 80)], "{order} {h:?}");
            } else {
                assert_eq!(atac_rad(&hits, false), [(200, 80)], "{order} {h:?}");
                assert_eq!(atac_rad(&hits, true), [(204, 71)], "{order} {h:?}");
            }
        }
    }
}

/// Contained mates mapped separately. (a) A 145-base reverse mate whose 3'
/// end runs 15 bases past a 100-base forward mate's 5' end: [200, 330), not
/// [185, 315). (b) A 120-base forward mate running past a 70-base reverse
/// mate's 5' end: [200, 300).
#[test]
fn contained_pair_spans_the_two_five_prime_ends() {
    let fx = fixture();
    let t = &fx.txp;
    let cases = [
        (
            "rc reads past fw",
            t[200..300].to_vec(),
            revcomp(&t[185..330]),
            200,
            130,
        ),
        (
            "fw reads past rc",
            t[200..320].to_vec(),
            revcomp(&t[230..300]),
            200,
            100,
        ),
    ];
    for (name, fw, rc, start, len) in cases {
        for (order, r1, r2) in [("fw/rc", &fw, &rc), ("rc/fw", &rc, &fw)] {
            let (mt, hits) = map_bulk(&fx.index, r1, r2);
            assert_eq!(mt, MappingType::MappedPair, "{name} {order}");
            assert_eq!(bulk_rad(&hits), [(start, len)], "{name} {order} {hits:?}");
            let (mt, hits) = map_atac(&fx.index, r1, r2);
            assert_eq!(mt, MappingType::MappedPair, "{name} {order}");
            assert_eq!(
                atac_rad(&hits, false),
                [(start, len)],
                "{name} {order} {hits:?}"
            );
        }
    }
}

/// A dovetailed pair at the transcript start: insert [0, 80) with 20 adapter
/// bases on each mate. The forward mate is at 0, so nothing is clipped, and
/// the reverse mate's projected start (-20) does not pull the start below 0.
/// With a forward mate hanging 10 bases off the start instead (pos = -10),
/// the clip shortens the fragment by those 10 bases, in both read orders.
#[test]
fn overhang_clip_uses_the_forward_mate() {
    let fx = fixture();
    let t = &fx.txp;
    let mut fw = t[0..80].to_vec();
    fw.extend(adapter(20, 3));
    let mut rc = revcomp(&t[0..80]);
    rc.extend(adapter(20, 4));
    for (order, r1, r2) in [("fw/rc", &fw, &rc), ("rc/fw", &rc, &fw)] {
        let (_, hits) = map_bulk(&fx.index, r1, r2);
        assert_eq!(bulk_rad(&hits), [(0, 80)], "{order} {hits:?}");
    }
    let mut fw = adapter(10, 5);
    fw.extend_from_slice(&t[0..70]);
    let rc = revcomp(&t[300..380]);
    for (order, r1, r2) in [("fw/rc", &fw, &rc), ("rc/fw", &rc, &fw)] {
        let (_, hits) = map_bulk(&fx.index, r1, r2);
        assert_eq!(hits[0].pair_start(), -10, "{order}");
        assert_eq!(bulk_rad(&hits), [(0, 380)], "{order} {hits:?}");
        let (_, hits) = map_atac(&fx.index, r1, r2);
        assert_eq!(atac_rad(&hits, false), [(0, 380)], "{order} {hits:?}");
    }
}

/// scATAC mates merged by `find_overlap` when the shorter read lies inside the
/// longer one. The merged fragment runs from the forward mate's 5' end to the
/// reverse mate's 5' end and is not cut to the shorter read: a 100-base read
/// inside a 130-base read over [200, 330) is [200, 330), and with the longer
/// read running 10 bases past the shorter read's 5' end (a 140-base read over
/// [190, 330)) it is still [200, 330). Both strands, both read orders.
#[test]
fn contained_overlap_merge_keeps_the_whole_fragment() {
    let fx = fixture();
    let t = &fx.txp;
    let cases = [
        // (shorter forward, longer reverse): fragment on the forward strand
        ("short fw", t[200..300].to_vec(), revcomp(&t[200..330])),
        (
            "short fw, long reads through",
            t[200..300].to_vec(),
            revcomp(&t[190..330]),
        ),
        // (shorter reverse, longer forward): fragment on the reverse strand
        ("short rc", revcomp(&t[230..330]), t[200..330].to_vec()),
        ("short rc, long reads through", revcomp(&t[230..330]), {
            let mut v = t[200..330].to_vec();
            v.extend(adapter(10, 6));
            v
        }),
    ];
    for (name, short, long) in cases {
        for (order, r1, r2) in [("short/long", &short, &long), ("long/short", &long, &short)] {
            let (ov, hits) = map_atac_overlap(&fx.index, r1, r2);
            assert_ne!(ov, OverlapType::NoOverlap, "{name} {order}");
            assert_eq!(hits.len(), 1, "{name} {order} {hits:?}");
            assert_eq!(hits[0].fragment_length, 130, "{name} {order} {hits:?}");
            assert_eq!(
                atac_rad(&hits, false),
                [(200, 130)],
                "{name} {order} {hits:?}"
            );
            assert_eq!(atac_rad(&hits, true), [(204, 121)], "{name} {order}");
        }
    }
}
