use std::path::Path;
fn main() {
    let cpp_rad = Path::new("/tmp/cpp_atac_out/map.rad");
    let rust_rad = Path::new("/tmp/rust_atac_out/map.rad");

    let result = piscem_rs::verify::rad_compare::compare_atac_rad_full(cpp_rad, rust_rad)
        .expect("failed to compare");

    println!("Header match: {}", result.header_match);
    println!("Records A (C++): {}", result.total_records_a);
    println!("Records B (Rust): {}", result.total_records_b);
    println!("Matching: {}", result.matching_records);
    println!("Missing in A: {}", result.missing_in_a);
    println!("Missing in B: {}", result.missing_in_b);
    println!("Notes: {}", result.notes);

    let match_rate = if result.total_records_a > 0 {
        result.matching_records as f64 / result.total_records_a as f64 * 100.0
    } else { 0.0 };
    println!("Match rate: {:.2}%", match_rate);

    if !result.first_mismatches.is_empty() {
        println!("First mismatches:");
        for m in &result.first_mismatches {
            println!("  - {}", m);
        }
    }
    println!("Same targets, diff detail: {}", result.same_targets_diff_detail);
    println!("Different targets: {}", result.different_targets);
    println!("Diff position only: {}", result.diff_pos_only);
    println!("Diff frag_len only: {}", result.diff_frag_len_only);
    println!("Diff num alignments: {}", result.diff_num_alns);
}
