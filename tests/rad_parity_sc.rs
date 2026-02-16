//! Single-cell RAD parity integration tests.
//!
//! These tests generate synthetic Chromium V3 reads from existing bulk
//! test data, run both the C++ and Rust SC mappers, and compare their
//! RAD output at the record level using multiset comparison.
//!
//! Requires test_data/ to be present. Run with:
//!   cargo test --features parity-test --release --test rad_parity_sc -- --ignored --nocapture

#![cfg(feature = "parity-test")]

use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

use anyhow::{Context, Result};

use piscem_rs::index::build::BuildConfig;
use piscem_rs::verify::rad_compare::compare_sc_rad_full;

// ---------------------------------------------------------------------------
// Paths
// ---------------------------------------------------------------------------

const CFISH_PREFIX: &str =
    "test_data/gencode_pc_v44_dbg/gencode_pc_v44_index_cfish";
const CPP_INDEX_PREFIX: &str =
    "test_data/gencode_pc_v44_index_cpp_from_dbg/gencode_pc_v44_index";
const RUST_INDEX_DIR: &str = "test_data/gencode_pc_v44_index_rust";
const RUST_INDEX_PREFIX: &str =
    "test_data/gencode_pc_v44_index_rust/gencode_pc_v44_index";
const BULK_READ2: &str = "test_data/sim_1M_2.fq.gz";
const CPP_SC_BIN: &str = "piscem-cpp/build/pesc-sc";
const CPP_BUILD_BIN: &str = "piscem-cpp/build/build";
const CPP_INDEX_DIR: &str = "test_data/gencode_pc_v44_index_cpp_from_dbg";
const SC_TEST_DIR: &str = "test_data/sc_test_data";

// Chromium V3: 16bp BC + 12bp UMI
const BC_LEN: usize = 16;
const UMI_LEN: usize = 12;
const R1_SUFFIX_LEN: usize = 22; // padding after BC+UMI to make R1 ~50bp

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Build the C++ index from cuttlefish output if it doesn't already exist.
fn ensure_cpp_index() -> Result<()> {
    let marker = PathBuf::from(CPP_INDEX_PREFIX).with_extension("sshash");
    if marker.exists() {
        eprintln!("C++ index already exists at {}", marker.display());
        return Ok(());
    }

    eprintln!("Building C++ index from cuttlefish output...");
    std::fs::create_dir_all(CPP_INDEX_DIR)
        .context("creating C++ index directory")?;

    let status = Command::new(CPP_BUILD_BIN)
        .arg("-i").arg(CFISH_PREFIX)
        .arg("-k").arg("31")
        .arg("-m").arg("19")
        .arg("-o").arg(CPP_INDEX_PREFIX)
        .arg("--canonical")
        .arg("--build-ec-table")
        .arg("-t").arg("4")
        .arg("--quiet")
        .status()
        .context("failed to run C++ index builder")?;
    if !status.success() {
        anyhow::bail!("C++ index builder failed with status: {}", status);
    }
    Ok(())
}

/// Build the Rust-format index if it doesn't already exist.
fn ensure_rust_index() -> Result<()> {
    let marker = PathBuf::from(RUST_INDEX_PREFIX).with_extension("ssi");
    if marker.exists() {
        eprintln!("Rust index already exists at {}", marker.display());
        return Ok(());
    }

    eprintln!("Building Rust index from cuttlefish output...");
    std::fs::create_dir_all(RUST_INDEX_DIR)
        .context("creating rust index directory")?;

    let config = BuildConfig {
        input_prefix: PathBuf::from(CFISH_PREFIX),
        output_prefix: PathBuf::from(RUST_INDEX_PREFIX),
        k: 31,
        m: 19,
        build_ec_table: true,
        num_threads: 0,
        canonical: true,
        seed: 1,
        single_mphf: false,
    };

    piscem_rs::index::build::build_index(&config)?;
    Ok(())
}

/// Generate a deterministic 16bp barcode from a seed index.
fn make_barcode(idx: usize) -> [u8; BC_LEN] {
    const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
    let mut bc = [b'A'; BC_LEN];
    let mut val = idx;
    for b in bc.iter_mut() {
        *b = BASES[val & 3];
        val >>= 2;
    }
    bc
}

/// Generate a deterministic 12bp UMI from a seed index.
fn make_umi(idx: usize) -> [u8; UMI_LEN] {
    const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
    let mut umi = [b'A'; UMI_LEN];
    // Use a simple hash to get different UMI per read
    let mut val = idx.wrapping_mul(2654435761); // Knuth multiplicative hash
    for b in umi.iter_mut() {
        *b = BASES[val & 3];
        val >>= 2;
    }
    umi
}

/// Generate synthetic Chromium V3 test data from existing bulk R2 reads.
///
/// R1 = [16bp BC][12bp UMI][suffix of Ts]
/// R2 = biological read (same as bulk)
///
/// Uses 100 rotating barcodes, deterministic UMIs.
fn generate_sc_test_data(num_reads: usize) -> Result<(PathBuf, PathBuf)> {
    std::fs::create_dir_all(SC_TEST_DIR)?;
    let r1_path = PathBuf::from(SC_TEST_DIR).join(format!("sc_v3_{}k_r1.fq.gz", num_reads / 1000));
    let r2_path = PathBuf::from(SC_TEST_DIR).join(format!("sc_v3_{}k_r2.fq.gz", num_reads / 1000));

    // Check if already cached
    if r1_path.exists() && r2_path.exists() {
        eprintln!("SC test data already exists at {}", SC_TEST_DIR);
        return Ok((r1_path, r2_path));
    }

    eprintln!("Generating synthetic V3 SC test data ({num_reads} reads)...");

    // Read R2 from bulk test data (use first `num_reads` reads)
    let (r2_reader, _) = niffler::send::from_path(BULK_READ2)
        .with_context(|| format!("open {}", BULK_READ2))?;
    let r2_bufreader = std::io::BufReader::new(r2_reader);

    // Parse R2 FASTQ records
    let mut r2_records: Vec<(Vec<u8>, Vec<u8>, Vec<u8>)> = Vec::new(); // (name, seq, qual)
    let mut lines = std::io::BufRead::lines(r2_bufreader);

    while r2_records.len() < num_reads {
        let name = match lines.next() {
            Some(Ok(l)) => l,
            _ => break,
        };
        let seq = match lines.next() {
            Some(Ok(l)) => l,
            _ => break,
        };
        let _plus = lines.next(); // skip +
        let qual = match lines.next() {
            Some(Ok(l)) => l,
            _ => break,
        };
        r2_records.push((name.into_bytes(), seq.into_bytes(), qual.into_bytes()));
    }

    if r2_records.len() < num_reads {
        anyhow::bail!(
            "Only {} reads in {}, need {}",
            r2_records.len(), BULK_READ2, num_reads
        );
    }

    // Write R1 (synthetic) and R2 (pass-through) as gzipped FASTQ
    let r1_file = std::fs::File::create(&r1_path)?;
    let r1_gz = flate2::write::GzEncoder::new(r1_file, flate2::Compression::fast());
    let mut r1_writer = std::io::BufWriter::new(r1_gz);

    let r2_file = std::fs::File::create(&r2_path)?;
    let r2_gz = flate2::write::GzEncoder::new(r2_file, flate2::Compression::fast());
    let mut r2_writer = std::io::BufWriter::new(r2_gz);

    let num_barcodes = 100;
    let suffix = vec![b'T'; R1_SUFFIX_LEN];
    let r1_qual = vec![b'I'; BC_LEN + UMI_LEN + R1_SUFFIX_LEN];

    for (i, (name, seq, qual)) in r2_records.iter().enumerate().take(num_reads) {
        let bc = make_barcode(i % num_barcodes);
        let umi = make_umi(i);

        // R1: @name\n[BC][UMI][suffix]\n+\n[qual]\n
        r1_writer.write_all(&name)?;
        r1_writer.write_all(b"\n")?;
        r1_writer.write_all(&bc)?;
        r1_writer.write_all(&umi)?;
        r1_writer.write_all(&suffix)?;
        r1_writer.write_all(b"\n+\n")?;
        r1_writer.write_all(&r1_qual)?;
        r1_writer.write_all(b"\n")?;

        // R2: pass-through
        r2_writer.write_all(&name)?;
        r2_writer.write_all(b"\n")?;
        r2_writer.write_all(seq)?;
        r2_writer.write_all(b"\n+\n")?;
        r2_writer.write_all(qual)?;
        r2_writer.write_all(b"\n")?;
    }

    r1_writer.flush()?;
    r2_writer.flush()?;
    drop(r1_writer);
    drop(r2_writer);

    eprintln!("SC test data generated: {} reads", num_reads);
    Ok((r1_path, r2_path))
}

/// Run the C++ SC mapper.
fn run_cpp_sc(
    r1: &Path,
    r2: &Path,
    output_dir: &Path,
    threads: usize,
) -> Result<()> {
    std::fs::create_dir_all(output_dir)?;
    let mut cmd = Command::new(CPP_SC_BIN);
    cmd.arg("-i").arg(CPP_INDEX_PREFIX)
        .arg("-o").arg(output_dir)
        .arg("-g").arg("chromium_v3")
        .arg("-t").arg(threads.to_string())
        .arg("-1").arg(r1)
        .arg("-2").arg(r2)
        .arg("--no-poison")
        .arg("--quiet");
    eprintln!("Running C++ SC mapper: {:?}", cmd);
    let status = cmd.status().context("failed to run C++ SC mapper")?;
    if !status.success() {
        anyhow::bail!("C++ SC mapper exited with status: {}", status);
    }
    Ok(())
}

/// Run the Rust SC mapper as a subprocess.
fn run_rust_sc(
    r1: &Path,
    r2: &Path,
    output_dir: &Path,
    threads: usize,
) -> Result<()> {
    let bin = PathBuf::from("target/release/piscem-rs");
    let bin = if bin.exists() {
        bin
    } else {
        PathBuf::from("target/debug/piscem-rs")
    };

    let mut cmd = Command::new(&bin);
    cmd.arg("map-scrna")
        .arg("-i").arg(RUST_INDEX_PREFIX)
        .arg("-o").arg(output_dir)
        .arg("-g").arg("chromium_v3")
        .arg("-t").arg(threads.to_string())
        .arg("-1").arg(r1)
        .arg("-2").arg(r2)
        .arg("--no-poison");
    eprintln!("Running Rust SC mapper: {:?}", cmd);
    let status = cmd.status().context("failed to run Rust SC mapper")?;
    if !status.success() {
        anyhow::bail!("Rust SC mapper exited with status: {}", status);
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

/// SC V3 RAD parity test: C++ vs Rust (10K reads, fast).
///
/// Run with: `cargo test --features parity-test --release --test rad_parity_sc -- --ignored --nocapture`
#[test]
#[ignore]
fn sc_v3_rad_parity() {
    // Check prerequisites
    if !Path::new(BULK_READ2).exists() {
        eprintln!("SKIP: {} not found", BULK_READ2);
        return;
    }
    if !Path::new(CPP_SC_BIN).exists() {
        eprintln!("SKIP: C++ SC binary not found at {}", CPP_SC_BIN);
        return;
    }
    if !Path::new(&format!("{CFISH_PREFIX}.cf_seg")).exists() {
        eprintln!("SKIP: cuttlefish output not found at {}", CFISH_PREFIX);
        return;
    }

    // 1. Build both indices
    ensure_cpp_index().expect("failed to build C++ index");
    ensure_rust_index().expect("failed to build Rust index");

    // 2. Generate synthetic SC test data (10K reads)
    let (r1_path, r2_path) = generate_sc_test_data(10_000)
        .expect("failed to generate SC test data");

    // 3. Create temp directories
    let tmpdir = tempfile::tempdir().expect("failed to create tempdir");
    let cpp_out_dir = tmpdir.path().join("cpp_out");
    let rust_out_dir = tmpdir.path().join("rust_out");

    // 4. Run both mappers with 1 thread for determinism
    run_cpp_sc(&r1_path, &r2_path, &cpp_out_dir, 1)
        .expect("C++ SC mapper failed");
    run_rust_sc(&r1_path, &r2_path, &rust_out_dir, 1)
        .expect("Rust SC mapper failed");

    // 5. Compare RAD files
    let cpp_rad = cpp_out_dir.join("map.rad");
    let rust_rad = rust_out_dir.join("map.rad");

    assert!(cpp_rad.exists(), "C++ RAD file not found: {}", cpp_rad.display());
    assert!(rust_rad.exists(), "Rust RAD file not found: {}", rust_rad.display());

    eprintln!("Comparing SC RAD files...");
    let result = compare_sc_rad_full(&cpp_rad, &rust_rad)
        .expect("failed to compare SC RAD files");

    eprintln!("SC Comparison result:");
    eprintln!("  Header match: {}", result.header_match);
    eprintln!("  Records A (C++): {}", result.total_records_a);
    eprintln!("  Records B (Rust): {}", result.total_records_b);
    eprintln!("  Matching: {}", result.matching_records);
    eprintln!("  Missing in A: {}", result.missing_in_a);
    eprintln!("  Missing in B: {}", result.missing_in_b);
    eprintln!("  Notes: {}", result.notes);
    if !result.first_mismatches.is_empty() {
        eprintln!("  First mismatches:");
        for m in &result.first_mismatches {
            eprintln!("    - {}", m);
        }
    }

    let match_rate = if result.total_records_a > 0 {
        result.matching_records as f64 / result.total_records_a as f64 * 100.0
    } else {
        0.0
    };
    eprintln!(
        "  C++ mapping rate: {:.2}% ({}/10K)",
        result.total_records_a as f64 / 100.0, result.total_records_a
    );
    eprintln!(
        "  Rust mapping rate: {:.2}% ({}/10K)",
        result.total_records_b as f64 / 100.0, result.total_records_b
    );
    eprintln!(
        "  Record match rate: {:.2}% ({}/{})",
        match_rate, result.matching_records, result.total_records_a
    );

    assert!(
        result.header_match,
        "RAD headers differ between C++ and Rust: {}",
        result.notes,
    );

    assert!(
        match_rate > 75.0,
        "Record-level match rate too low ({:.2}%): {}",
        match_rate,
        result.notes,
    );
}

/// SC V3 RAD parity test with 1M reads (slower, more thorough).
#[test]
#[ignore]
fn sc_v3_rad_parity_1m() {
    if !Path::new(BULK_READ2).exists() {
        eprintln!("SKIP: {} not found", BULK_READ2);
        return;
    }
    if !Path::new(CPP_SC_BIN).exists() {
        eprintln!("SKIP: C++ SC binary not found at {}", CPP_SC_BIN);
        return;
    }
    if !Path::new(&format!("{CFISH_PREFIX}.cf_seg")).exists() {
        eprintln!("SKIP: cuttlefish output not found at {}", CFISH_PREFIX);
        return;
    }

    ensure_cpp_index().expect("failed to build C++ index");
    ensure_rust_index().expect("failed to build Rust index");

    let (r1_path, r2_path) = generate_sc_test_data(1_000_000)
        .expect("failed to generate SC test data");

    let tmpdir = tempfile::tempdir().expect("failed to create tempdir");
    let cpp_out_dir = tmpdir.path().join("cpp_out");
    let rust_out_dir = tmpdir.path().join("rust_out");

    run_cpp_sc(&r1_path, &r2_path, &cpp_out_dir, 1)
        .expect("C++ SC mapper failed");
    run_rust_sc(&r1_path, &r2_path, &rust_out_dir, 1)
        .expect("Rust SC mapper failed");

    let cpp_rad = cpp_out_dir.join("map.rad");
    let rust_rad = rust_out_dir.join("map.rad");

    let result = compare_sc_rad_full(&cpp_rad, &rust_rad)
        .expect("failed to compare SC RAD files");

    eprintln!("SC 1M Comparison result:");
    eprintln!("  Records A (C++): {}", result.total_records_a);
    eprintln!("  Records B (Rust): {}", result.total_records_b);
    eprintln!("  Matching: {}", result.matching_records);
    eprintln!("  Notes: {}", result.notes);

    let match_rate = if result.total_records_a > 0 {
        result.matching_records as f64 / result.total_records_a as f64 * 100.0
    } else {
        0.0
    };
    eprintln!("  Record match rate: {:.2}%", match_rate);

    assert!(
        match_rate > 75.0,
        "Record-level match rate too low ({:.2}%)",
        match_rate,
    );
}
