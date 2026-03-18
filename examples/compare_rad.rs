/// Quick tool to semantically compare two RAD files.
///
/// Usage:
///   cargo run --example compare_rad --features parity-test -- bulk <file_a> <file_b>
///   cargo run --example compare_rad --features parity-test -- sc <file_a> <file_b>
///   cargo run --example compare_rad --features parity-test -- atac <file_a> <file_b>
///
/// Handles reordered records (multithreading) and reordered reference headers.
use std::path::Path;

use piscem_rs::verify::rad_compare;

fn usage() -> ! {
    eprintln!("Usage: compare_rad <bulk|sc|atac> <file_a> <file_b>");
    std::process::exit(1);
}

fn print_summary(result: &rad_compare::RecordComparisonSummary) {
    println!("Header match:      {}", result.header_match);
    println!("Records in A:      {}", result.total_records_a);
    println!("Records in B:      {}", result.total_records_b);
    println!("Matching:          {}", result.matching_records);
    println!("Missing in A:      {}", result.missing_in_a);
    println!("Missing in B:      {}", result.missing_in_b);
    println!("PASSED:            {}", result.passed);
    println!("Notes:             {}", result.notes);

    let denom = result.total_records_a.max(result.total_records_b);
    let match_rate = if denom > 0 {
        result.matching_records as f64 / denom as f64 * 100.0
    } else {
        100.0
    };
    println!("Match rate:        {match_rate:.4}%");

    if !result.first_mismatches.is_empty() {
        println!("\nFirst mismatches (up to 10):");
        for m in &result.first_mismatches {
            println!("  - {m}");
        }
    }

    // Bulk-only detail fields
    if result.same_targets_diff_detail > 0 || result.different_targets > 0 {
        println!("\nMismatch breakdown:");
        println!(
            "  Same targets, diff detail:  {}",
            result.same_targets_diff_detail
        );
        println!("    Diff position only:       {}", result.diff_pos_only);
        println!(
            "    Diff frag_len only:       {}",
            result.diff_frag_len_only
        );
        println!("    Diff num alignments:      {}", result.diff_num_alns);
        println!("  Different targets entirely: {}", result.different_targets);
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 4 {
        usage();
    }

    let mode = args[1].as_str();
    let path_a = Path::new(&args[2]);
    let path_b = Path::new(&args[3]);

    println!("Comparing RAD files ({mode}):");
    println!("  A: {}", path_a.display());
    println!("  B: {}", path_b.display());
    println!();

    let result = match mode {
        "bulk" => rad_compare::compare_bulk_rad_full(path_a, path_b),
        "sc" => rad_compare::compare_sc_rad_full(path_a, path_b),
        "atac" => rad_compare::compare_atac_rad_full(path_a, path_b),
        _ => {
            eprintln!("Unknown mode '{mode}'. Must be bulk, sc, or atac.");
            usage();
        }
    };

    match result {
        Ok(summary) => {
            print_summary(&summary);
            if !summary.passed {
                std::process::exit(1);
            }
        }
        Err(e) => {
            eprintln!("Error: {e:#}");
            std::process::exit(2);
        }
    }
}
