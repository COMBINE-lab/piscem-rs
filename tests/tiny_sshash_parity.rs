//! Tiny vs. sshash parity tests.
//!
//! Loads a small Rust-format index, builds a parallel Tiny-backed view, and
//! asserts that both dictionaries + contig tables return identical
//! mapping-relevant information. Skipped if `test_data/gencode_pc_v44_index_rust`
//! is not present.
//!
//! Run with: `cargo test --test tiny_sshash_parity -- --nocapture`

use piscem_rs::index::contig_table::{ContigTableLike, TinyContigTable};
use piscem_rs::index::reference_index::ReferenceIndex;
use sshash_lib::{KmerDictionary, KmerStreamingQuery, LookupResult, dispatch_on_k};
use tiny_dict::TinyDictionary;

/// Default paths — override with `PISCEM_PARITY_INDEX` and
/// `PISCEM_PARITY_FASTQ` env vars. If neither file exists the test is skipped.
const DEFAULT_INDEX_PREFIX: &str = "test_data/gencode_pc_v44_index_rust/gencode_pc_v44_index";
const DEFAULT_TEST_FASTQ: &str = "test_data/sim_1M_1.fq.gz";
const MAX_READS: usize = 2000;

fn index_prefix() -> std::path::PathBuf {
    std::env::var("PISCEM_PARITY_INDEX")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|_| std::path::PathBuf::from(DEFAULT_INDEX_PREFIX))
}

fn test_fastq() -> std::path::PathBuf {
    std::env::var("PISCEM_PARITY_FASTQ")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|_| std::path::PathBuf::from(DEFAULT_TEST_FASTQ))
}

/// Fields compared across dictionary implementations. `kmer_id` and
/// `kmer_offset` are storage-order-dependent (different representations may
/// number k-mers differently) and are intentionally excluded — the mapping
/// pipeline consumes `string_id` / `kmer_id_in_string` / `string_begin` /
/// `string_end` / `kmer_orientation`, which must agree.
#[derive(Debug, PartialEq, Eq)]
struct CompareFields {
    is_found: bool,
    kmer_orientation: i64,
    kmer_id_in_string: u64,
    string_id: u64,
    string_begin: u64,
    string_end: u64,
}

impl From<&LookupResult> for CompareFields {
    fn from(r: &LookupResult) -> Self {
        CompareFields {
            is_found: r.is_found(),
            kmer_orientation: if r.is_found() { r.kmer_orientation } else { 0 },
            kmer_id_in_string: if r.is_found() { r.kmer_id_in_string } else { 0 },
            string_id: if r.is_found() { r.string_id } else { 0 },
            string_begin: if r.is_found() { r.string_begin } else { 0 },
            string_end: if r.is_found() { r.string_end } else { 0 },
        }
    }
}

#[test]
fn tiny_sshash_lookup_parity_on_test_reads() {
    let index_prefix = index_prefix();
    if !index_prefix.with_extension("ssi").exists() {
        eprintln!("Skipping: sshash index not found at {}", index_prefix.display());
        return;
    }
    let index_prefix = index_prefix.as_path();
    let read1 = test_fastq();
    if !read1.exists() {
        eprintln!("Skipping: {} not found", read1.display());
        return;
    }
    let read1 = read1.as_path();

    let index = ReferenceIndex::load(index_prefix, false, false).expect("load sshash index");
    let k = index.k();
    eprintln!(
        "Loaded sshash index: k={}, {} strings, {} contigs",
        k,
        index.dict().num_strings(),
        index.contig_table().num_contigs()
    );

    dispatch_on_k!(k, K => {
        let tiny_dict = TinyDictionary::from_sshash::<K>(index.dict());

        let mut sshash_q = index.dict().create_streaming_query::<K>();
        let mut tiny_q = tiny_dict.create_streaming_query::<K>();

        let (reader, _) = niffler::send::from_path(read1).expect("open test fastq");
        let mut rdr = paraseq::fastq::Reader::new(reader);
        let mut rset = rdr.new_record_set();
        rset.fill(&mut rdr).expect("read record set");

        let mut compared = 0usize;
        let mut found = 0usize;
        use paraseq::Record;
        for (i, rec) in rset.iter().enumerate() {
            if i >= MAX_READS { break; }
            let rec = rec.expect("decode record");
            let seq = rec.seq();
            if seq.len() < K { continue; }

            sshash_q.reset();
            tiny_q.reset();

            for start in 0..=(seq.len() - K) {
                let kb = &seq[start..start+K];
                if kb.contains(&b'N') {
                    sshash_q.reset();
                    tiny_q.reset();
                    continue;
                }
                let a = sshash_q.lookup(kb);
                let b = tiny_q.lookup(kb);
                let af = CompareFields::from(&a);
                let bf = CompareFields::from(&b);
                assert_eq!(
                    af, bf,
                    "LookupResult mismatch on read {} pos {}: sshash={:?} tiny={:?}",
                    i, start, af, bf
                );
                compared += 1;
                if a.is_found() { found += 1; }
            }
        }
        eprintln!("compared {} k-mers, {} found ({:.1}%)",
            compared, found,
            if compared > 0 { 100.0 * found as f64 / compared as f64 } else { 0.0 });
        assert!(compared > 0, "no k-mers compared — check test data");
        assert!(found > 0, "no k-mers found — dictionary parity path is broken");
    });
}

#[test]
fn tiny_sshash_contig_entries_parity() {
    let index_prefix = index_prefix();
    if !index_prefix.with_extension("ssi").exists() {
        eprintln!("Skipping: sshash index not found at {}", index_prefix.display());
        return;
    }
    let index_prefix = index_prefix.as_path();

    let index = ReferenceIndex::load(index_prefix, false, false).expect("load sshash index");
    let tiny_ct = TinyContigTable::from_contig_table(index.contig_table());

    let n = index.contig_table().num_contigs();
    assert_eq!(n, tiny_ct.num_contigs(), "contig count mismatch");

    let mut total_entries = 0usize;
    for cid in 0..(n as u64) {
        let a = index.contig_table().contig_entries(cid);
        let b = tiny_ct.contig_entries(cid);
        assert_eq!(a.len(), b.len(), "entry count mismatch for contig {}", cid);

        // Compare as multisets: entry ordering within a contig is not part of
        // the observable contract, but the decoded (ref_id, pos, orient) tuples
        // must be identical. Entries from `contig_entries` are already the
        // packed u64 encoding under the same EntryEncoding, so direct value
        // comparison is sufficient.
        let mut av: Vec<u64> = a.iter().collect();
        let mut bv: Vec<u64> = b.iter().collect();
        av.sort_unstable();
        bv.sort_unstable();
        assert_eq!(av, bv, "entry set mismatch for contig {}", cid);
        total_entries += a.len();
    }
    eprintln!(
        "compared {} contigs, {} total entries; tiny inline fraction = {:.2}%",
        n,
        total_entries,
        tiny_ct.inline_fraction() * 100.0
    );
}
