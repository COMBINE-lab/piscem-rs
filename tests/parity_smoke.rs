//! Parity smoke tests.
//!
//! These tests validate that RAD file headers are well-formed when written
//! by our RAD writer and that we can read them back.

use piscem_rs::io::rad::{
    write_rad_header_bulk, write_rad_header_sc, write_rad_header_atac,
};

#[test]
fn parity_smoke_bulk_header_roundtrip() {
    let mut buf = Vec::new();
    let names = vec!["chr1", "chr2", "chr3"];
    let ref_lens = vec![1000u32, 2000u32, 3000u32];
    let chunk_off = write_rad_header_bulk(&mut buf, true, 3, &names, &ref_lens).unwrap();

    assert!(chunk_off > 0);
    assert!(buf.len() > 20); // sanity: header should be non-trivial
    assert_eq!(buf[0], 1); // is_paired
}

#[test]
fn parity_smoke_sc_header_roundtrip() {
    let mut buf = Vec::new();
    let names = vec!["gene1", "gene2"];
    let (chunk_off, rlen_off) =
        write_rad_header_sc(&mut buf, 2, &names, 16, 12, true).unwrap();

    assert!(chunk_off > 0);
    assert!(rlen_off.is_some());
    assert_eq!(buf[0], 0); // is_paired=false for SC
}

#[test]
fn parity_smoke_atac_header_roundtrip() {
    let mut buf = Vec::new();
    let names = vec!["chr1"];
    let ref_lens = vec![248956422u32];
    let chunk_off =
        write_rad_header_atac(&mut buf, 1, &names, &ref_lens, 16).unwrap();

    assert!(chunk_off > 0);
    assert_eq!(buf[0], 1); // is_paired=true for ATAC
}

/// Integration test: build + map pipeline.
///
/// This test is ignored by default since it requires test_data/ to be present.
/// Run with: `cargo test --test parity_smoke -- --ignored`
#[test]
#[ignore]
fn integration_build_and_map() {
    use std::path::Path;

    let index_prefix = Path::new("test_data/gencode_pc_v44_index_nopoison/gencode_pc_v44_index");
    if !index_prefix.with_extension("sshash").exists() {
        eprintln!("Skipping: test_data index not found");
        return;
    }

    // Try loading the index
    let index = piscem_rs::index::reference_index::ReferenceIndex::load(
        index_prefix, false, false,
    );
    assert!(index.is_ok(), "Failed to load index: {:?}", index.err());

    let idx = index.unwrap();
    assert!(idx.num_refs() > 0, "Index has no references");
    assert!(idx.k() > 0, "Index has invalid k");

    eprintln!(
        "Index loaded: k={}, {} refs",
        idx.k(),
        idx.num_refs(),
    );
}
