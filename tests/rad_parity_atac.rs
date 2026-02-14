//! ATAC RAD parity integration tests.
//!
//! These tests run both the C++ and Rust scATAC mappers on the same reads
//! and compare their RAD output at the record level using multiset comparison.
//!
//! Requires test_data/atacseq/ to be present. Run with:
//!   cargo test --features parity-test --release --test rad_parity_atac -- --ignored --nocapture

// Only compile when the parity-test feature is enabled
#![cfg(feature = "parity-test")]

use std::path::{Path, PathBuf};
use std::process::Command;

use anyhow::{Context, Result};

use piscem_rs::verify::rad_compare::compare_atac_rad_full;

// ---------------------------------------------------------------------------
// Paths
// ---------------------------------------------------------------------------

const CPP_INDEX_PREFIX: &str = "test_data/atacseq/cpp_index/k25_m17";
const RUST_INDEX_PREFIX: &str = "test_data/atacseq/rust_index/k25_m17";

// Triple-file FASTQ: R1 (genomic left), R2 (barcode), R3 (genomic right)
const READ1: &str = "test_data/atacseq/atac_seq_5M_r1.fq.gz";
const BARCODE: &str = "test_data/atacseq/atac_seq_5M_r2.fq.gz";
const READ2: &str = "test_data/atacseq/atac_seq_5M_r3.fq.gz";

const CPP_ATAC_BIN: &str = "piscem-cpp/build/pesc-sc-atac";

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Run the C++ scATAC mapper.
fn run_cpp_atac(output_dir: &Path, threads: usize) -> Result<()> {
    // --no-poison defaults to 1 (true) in C++, so we don't pass it explicitly.
    // --quiet minimizes console output.
    let cmd = Command::new(CPP_ATAC_BIN)
        .arg("-i").arg(CPP_INDEX_PREFIX)
        .arg("-1").arg(READ1)
        .arg("-2").arg(READ2)
        .arg("-b").arg(BARCODE)
        .arg("-o").arg(output_dir)
        .arg("-t").arg(threads.to_string())
        .arg("--quiet")
        .status()
        .context("failed to run C++ ATAC mapper")?;
    if !cmd.success() {
        anyhow::bail!("C++ ATAC mapper exited with status: {}", cmd);
    }
    Ok(())
}

/// Run the Rust scATAC mapper as a subprocess.
fn run_rust_atac(output_dir: &Path, threads: usize) -> Result<()> {
    let bin = PathBuf::from("target/release/piscem-rs");
    let bin = if bin.exists() {
        bin
    } else {
        PathBuf::from("target/debug/piscem-rs")
    };

    let cmd = Command::new(&bin)
        .arg("map-scatac")
        .arg("-i").arg(RUST_INDEX_PREFIX)
        .arg("-1").arg(READ1)
        .arg("-2").arg(READ2)
        .arg("-b").arg(BARCODE)
        .arg("-o").arg(output_dir)
        .arg("-t").arg(threads.to_string())
        .status()
        .context("failed to run Rust ATAC mapper")?;
    if !cmd.success() {
        anyhow::bail!("Rust ATAC mapper exited with status: {}", cmd);
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

/// Full ATAC RAD parity test: C++ vs Rust.
///
/// Run with: `cargo test --features parity-test --release --test rad_parity_atac -- --ignored --nocapture`
#[test]
#[ignore]
fn atac_rad_parity() {
    // Check prerequisites
    if !Path::new(READ1).exists() {
        eprintln!("SKIP: {} not found", READ1);
        return;
    }
    if !Path::new(CPP_ATAC_BIN).exists() {
        eprintln!("SKIP: C++ binary not found at {}", CPP_ATAC_BIN);
        return;
    }
    if !Path::new(&format!("{CPP_INDEX_PREFIX}.sshash")).exists() {
        eprintln!("SKIP: C++ index not found at {}", CPP_INDEX_PREFIX);
        return;
    }
    if !Path::new(&format!("{RUST_INDEX_PREFIX}.ssi")).exists() {
        eprintln!("SKIP: Rust index not found at {}", RUST_INDEX_PREFIX);
        return;
    }

    // Create temp directories
    let tmpdir = tempfile::tempdir().expect("failed to create tempdir");
    let cpp_out_dir = tmpdir.path().join("cpp_out");
    let rust_out_dir = tmpdir.path().join("rust_out");
    std::fs::create_dir_all(&cpp_out_dir).unwrap();
    std::fs::create_dir_all(&rust_out_dir).unwrap();

    // Run both mappers with 1 thread for determinism
    eprintln!("Running C++ ATAC mapper...");
    run_cpp_atac(&cpp_out_dir, 1).expect("C++ ATAC mapper failed");

    eprintln!("Running Rust ATAC mapper...");
    run_rust_atac(&rust_out_dir, 1).expect("Rust ATAC mapper failed");

    // Compare RAD files
    let cpp_rad = cpp_out_dir.join("map.rad");
    let rust_rad = rust_out_dir.join("map.rad");

    assert!(cpp_rad.exists(), "C++ RAD file not found: {}", cpp_rad.display());
    assert!(rust_rad.exists(), "Rust RAD file not found: {}", rust_rad.display());

    eprintln!("Comparing ATAC RAD files...");
    let result = compare_atac_rad_full(&cpp_rad, &rust_rad)
        .expect("failed to compare ATAC RAD files");

    eprintln!("Comparison result:");
    eprintln!("  Header match: {}", result.header_match);
    eprintln!("  Records A (C++): {}", result.total_records_a);
    eprintln!("  Records B (Rust): {}", result.total_records_b);
    eprintln!("  Matching: {}", result.matching_records);
    eprintln!("  Missing in A: {}", result.missing_in_a);
    eprintln!("  Missing in B: {}", result.missing_in_b);
    eprintln!("  Mismatch categories:");
    eprintln!("    Same targets, diff detail: {}", result.same_targets_diff_detail);
    eprintln!("    Different targets: {}", result.different_targets);
    eprintln!("    Diff position only: {}", result.diff_pos_only);
    eprintln!("    Diff frag_len only: {}", result.diff_frag_len_only);
    eprintln!("    Diff num alignments: {}", result.diff_num_alns);
    eprintln!("  Notes: {}", result.notes);
    if !result.first_mismatches.is_empty() {
        eprintln!("  First mismatches:");
        for m in &result.first_mismatches {
            eprintln!("    - {}", m);
        }
    }

    // Report record-level parity
    let match_rate = if result.total_records_a > 0 {
        result.matching_records as f64 / result.total_records_a as f64 * 100.0
    } else {
        0.0
    };
    eprintln!("  Record match rate: {:.2}% ({}/{})",
        match_rate, result.matching_records, result.total_records_a);

    // Assert reasonable record-level parity
    assert!(
        match_rate > 50.0,
        "Record-level match rate too low ({:.2}%): {}",
        match_rate,
        result.notes,
    );
}
