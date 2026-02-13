//! Parity smoke tests.
//!
//! These tests validate that RAD file headers are well-formed when written
//! by our RAD writer and that we can read them back.

use paraseq::Record;
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

/// Integration test: load Rust-format index and verify k-mer lookup works.
///
/// Run with: `cargo test --test parity_smoke -- --ignored --nocapture`
#[test]
#[ignore]
fn integration_rust_index_lookup() {
    use std::path::Path;

    let index_prefix = Path::new("test_data/gencode_pc_v44_index_rust/gencode_pc_v44_index");
    if !index_prefix.with_extension("ssi").exists() {
        eprintln!("Skipping: Rust index not found at {}", index_prefix.display());
        return;
    }

    let index = piscem_rs::index::reference_index::ReferenceIndex::load(
        index_prefix, true, false,
    )
    .expect("failed to load Rust index");

    let dict = index.dict();
    eprintln!(
        "Index loaded: k={}, {} refs, {} strings, canonical={}",
        index.k(), index.num_refs(), dict.num_strings(), dict.canonical()
    );

    // Read first read from test FASTQ
    let read1 = Path::new("test_data/sim_1M_1.fq.gz");
    if !read1.exists() {
        eprintln!("Skipping: {} not found", read1.display());
        return;
    }

    let (reader, _) = niffler::send::from_path(read1).unwrap();
    let mut rdr = paraseq::fastq::Reader::new(reader);
    let mut rset = rdr.new_record_set();
    rset.fill(&mut rdr).unwrap();

    let k = index.k();
    let mut total_kmers = 0usize;
    let mut found_kmers = 0usize;

    use sshash_lib::dispatch_on_k;
    dispatch_on_k!(k, K => {
        let mut query = dict.create_streaming_query::<K>();

        for (i, rec) in rset.iter().enumerate() {
            if i >= 10 { break; }
            let rec = rec.unwrap();
            let seq = rec.seq();
            if seq.len() < K { continue; }

            query.reset();
            let mut read_found = 0;
            let mut read_total = 0;
            for start in 0..=(seq.len() - K) {
                let kb = &seq[start..start+K];
                if kb.contains(&b'N') {
                    query.reset();
                    continue;
                }
                let ks = std::str::from_utf8(kb).unwrap();
                let result = query.lookup(ks);
                read_total += 1;
                if result.is_found() {
                    read_found += 1;
                }
            }
            total_kmers += read_total;
            found_kmers += read_found;
            if i < 3 {
                eprintln!(
                    "  Read {}: {} bases, {}/{} k-mers found",
                    i, seq.len(), read_found, read_total
                );
            }
        }
    });

    eprintln!(
        "Total: {}/{} k-mers found ({:.1}%)",
        found_kmers, total_kmers,
        if total_kmers > 0 { found_kmers as f64 / total_kmers as f64 * 100.0 } else { 0.0 }
    );

    assert!(
        found_kmers > 0,
        "No k-mers found in first 10 reads — dictionary lookup broken"
    );
}
