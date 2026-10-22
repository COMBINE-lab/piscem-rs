//! Parity smoke tests.
//!
//! These tests validate that RAD file headers are well-formed when written
//! by our RAD writer and that libradicl (alevin-fry's reader) reads them back.

use libradicl::header::RadPrelude;
use paraseq::Record;
use piscem_rs::io::rad::{write_rad_header_atac, write_rad_header_bulk, write_rad_header_sc};

/// Parse a header we wrote back with libradicl's reader (the one alevin-fry
/// uses) and check the fields every RAD mode shares.
fn read_back(buf: &[u8], is_paired: u8, names: &[&str]) -> RadPrelude {
    let p = RadPrelude::from_bytes(&mut &buf[..]).expect("libradicl rejects our header");
    assert_eq!(&buf[..8], b"RAD_FILE", "RAD v2 magic");
    assert_eq!(p.hdr.version.major(), 2);
    assert_eq!(p.hdr.is_paired, is_paired);
    assert_eq!(p.hdr.ref_count, names.len() as u64);
    assert_eq!(p.hdr.ref_names, names);
    p
}

#[test]
fn parity_smoke_bulk_header_roundtrip() {
    let mut buf = Vec::new();
    let names = vec!["chr1", "chr2", "chr3"];
    let ref_lens = vec![1000u32, 2000u32, 3000u32];
    let chunk_off = write_rad_header_bulk(&mut buf, true, 3, &names, &ref_lens).unwrap();

    assert!(chunk_off > 0);
    read_back(&buf, 1, &names);
}

#[test]
fn parity_smoke_sc_header_roundtrip() {
    let mut buf = Vec::new();
    let names = vec!["gene1", "gene2"];
    let (chunk_off, rlen_off) = write_rad_header_sc(&mut buf, 2, &names, 16, 12, true).unwrap();

    assert!(chunk_off > 0);
    assert!(rlen_off.is_some());
    let p = read_back(&buf, 0, &names);
    assert_eq!(p.read_tags.tags.len(), 2); // barcode, UMI
}

#[test]
fn parity_smoke_atac_header_roundtrip() {
    let mut buf = Vec::new();
    let names = vec!["chr1"];
    let ref_lens = vec![248956422u32];
    let chunk_off = write_rad_header_atac(&mut buf, 1, &names, &ref_lens, 16).unwrap();

    assert!(chunk_off > 0);
    read_back(&buf, 1, &names);
}

/// Integration test: load Rust-format index and verify k-mer lookup works.
///
/// Run with: `cargo test --test parity_smoke -- --ignored --nocapture`
#[test]
#[ignore]
fn integration_rust_index_lookup() {
    use piscem_rs::index::paths::{DICT_EXT, index_artifact_path};
    use std::path::Path;

    let index_prefix = Path::new("test_data/gencode_pc_v44_index_rust/gencode_pc_v44_index");
    if !index_artifact_path(index_prefix, DICT_EXT).exists() {
        eprintln!(
            "Skipping: Rust index not found at {}",
            index_prefix.display()
        );
        return;
    }

    let index = piscem_rs::index::reference_index::ReferenceIndex::load(index_prefix, true, false)
        .expect("failed to load Rust index");

    let dict = index.dict();
    eprintln!(
        "Index loaded: k={}, {} refs, {} strings",
        index.k(),
        index.num_refs(),
        dict.num_strings()
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
                let result = query.lookup(kb);
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
        found_kmers,
        total_kmers,
        if total_kmers > 0 {
            found_kmers as f64 / total_kmers as f64 * 100.0
        } else {
            0.0
        }
    );

    assert!(
        found_kmers > 0,
        "No k-mers found in first 10 reads — dictionary lookup broken"
    );
}
