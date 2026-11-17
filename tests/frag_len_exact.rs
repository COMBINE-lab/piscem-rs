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
use piscem_rs::mapping::map_fragment::{map_pe_fragment, map_pe_fragment_atac};
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
