/// Dump basic index statistics for debugging.
use std::path::Path;
use piscem_rs::index::reference_index::ReferenceIndex;

fn main() {
    let prefix = std::env::args().nth(1).expect("usage: index_stats <index_prefix>");
    let index = ReferenceIndex::load(Path::new(&prefix), false, false)
        .expect("failed to load index");

    println!("k = {}", index.k());
    println!("num_refs = {}", index.num_refs());
    let total_len: u64 = (0..index.num_refs()).map(|i| index.ref_len(i)).sum();
    println!("total_ref_length = {}", total_len);

    // Print first 10 and last 5 ref lengths
    let n = index.num_refs();
    println!("\nFirst 10 references:");
    for i in 0..10.min(n) {
        println!("  ref[{}] = {} (len={})", i, index.ref_name(i), index.ref_len(i));
    }
    if n > 15 {
        println!("  ...");
        println!("Last 5 references:");
        for i in (n-5)..n {
            println!("  ref[{}] = {} (len={})", i, index.ref_name(i), index.ref_len(i));
        }
    }

    // Contig table stats
    let ctab = index.contig_table();
    println!("\nContig table:");
    println!("  num_contigs = {}", ctab.num_contigs());
    println!("  num_entries = {}", ctab.num_entries());
    // Sample a few k-mer lookups
    let dict = index.dict();
    println!("  dict num_strings = {}", dict.num_strings());
    println!("  dict k = {}", dict.k());

    // Sample contig entry counts to compare index density
    let mut total_entries_check: usize = 0;
    let mut max_entries_per_contig: usize = 0;
    let mut contigs_with_zero_entries: usize = 0;
    for c in 0..ctab.num_contigs() {
        let entries = ctab.contig_entries(c as u64);
        let count = entries.len();
        total_entries_check += count;
        if count > max_entries_per_contig {
            max_entries_per_contig = count;
        }
        if count == 0 {
            contigs_with_zero_entries += 1;
        }
    }
    println!("  total_entries (verified) = {}", total_entries_check);
    println!("  max_entries_per_contig = {}", max_entries_per_contig);
    println!("  contigs_with_zero_entries = {}", contigs_with_zero_entries);
}
