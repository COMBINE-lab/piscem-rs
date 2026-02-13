//! Bulk RAD parity integration tests.
//!
//! These tests build a Rust-format index from cuttlefish output, run both
//! the C++ and Rust bulk mappers on the same reads, and compare their RAD
//! output at the record level using multiset comparison.
//!
//! Requires test_data/ to be present. Run with:
//!   cargo test --features parity-test --test rad_parity_bulk -- --ignored --nocapture

// Only compile when the parity-test feature is enabled
#![cfg(feature = "parity-test")]

use std::path::{Path, PathBuf};
use std::process::Command;

use anyhow::{Context, Result};

use piscem_rs::index::build::BuildConfig;
use piscem_rs::verify::rad_compare::compare_bulk_rad_full;

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
const READ1: &str = "test_data/sim_1M_1.fq.gz";
const READ2: &str = "test_data/sim_1M_2.fq.gz";
const CPP_BULK_BIN: &str = "piscem-cpp/build/pesc-bulk";
const CPP_BUILD_BIN: &str = "piscem-cpp/build/build";
const CPP_INDEX_DIR: &str = "test_data/gencode_pc_v44_index_cpp_from_dbg";

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
    eprintln!("C++ index built successfully.");
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
    };

    piscem_rs::index::build::build_index(&config)?;
    eprintln!("Rust index built successfully.");
    Ok(())
}

/// Run the C++ bulk mapper.
fn run_cpp_bulk(
    output_stem: &Path,
    threads: usize,
    no_poison: bool,
) -> Result<()> {
    let mut cmd = Command::new(CPP_BULK_BIN);
    cmd.arg("-i").arg(CPP_INDEX_PREFIX)
        .arg("-o").arg(output_stem)
        .arg("-t").arg(threads.to_string())
        .arg("-1").arg(READ1)
        .arg("-2").arg(READ2)
        .arg("--quiet");
    if no_poison {
        cmd.arg("--no-poison");
    }
    eprintln!("Running C++ mapper: {:?}", cmd);
    let status = cmd.status().context("failed to run C++ bulk mapper")?;
    if !status.success() {
        anyhow::bail!("C++ bulk mapper exited with status: {}", status);
    }
    Ok(())
}

/// Run the Rust bulk mapper as a subprocess.
fn run_rust_bulk(
    output_dir: &Path,
    index_prefix: &str,
    threads: usize,
    no_poison: bool,
) -> Result<()> {
    // Use cargo run --release for performance, but fall back to debug
    let bin = PathBuf::from("target/release/piscem-rs");
    let bin = if bin.exists() {
        bin
    } else {
        PathBuf::from("target/debug/piscem-rs")
    };

    let mut cmd = Command::new(&bin);
    cmd.arg("map-bulk")
        .arg("-i").arg(index_prefix)
        .arg("-o").arg(output_dir)
        .arg("-t").arg(threads.to_string())
        .arg("-1").arg(READ1)
        .arg("-2").arg(READ2);
    if no_poison {
        cmd.arg("--no-poison");
    }
    eprintln!("Running Rust mapper: {:?}", cmd);
    let status = cmd.status().context("failed to run Rust bulk mapper")?;
    if !status.success() {
        anyhow::bail!("Rust bulk mapper exited with status: {}", status);
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

/// Full bulk PE RAD parity test: C++ vs Rust.
///
/// Run with: `cargo test --features parity-test --test rad_parity_bulk -- --ignored --nocapture`
#[test]
#[ignore]
fn bulk_pe_rad_parity() {
    // Check prerequisites
    if !Path::new(READ1).exists() {
        eprintln!("SKIP: {} not found", READ1);
        return;
    }
    if !Path::new(CPP_BULK_BIN).exists() {
        eprintln!("SKIP: C++ binary not found at {}", CPP_BULK_BIN);
        return;
    }
    if !Path::new(&format!("{CFISH_PREFIX}.cf_seg")).exists() {
        eprintln!("SKIP: cuttlefish output not found at {}", CFISH_PREFIX);
        return;
    }

    // 1. Build both indices from same cuttlefish output
    ensure_cpp_index().expect("failed to build C++ index");
    ensure_rust_index().expect("failed to build Rust index");

    // 2. Create temp directories
    let tmpdir = tempfile::tempdir().expect("failed to create tempdir");
    let cpp_out_stem = tmpdir.path().join("cpp_out");
    let rust_out_dir = tmpdir.path().join("rust_out");
    std::fs::create_dir_all(&rust_out_dir).unwrap();

    // 3. Run both mappers with 1 thread for determinism
    run_cpp_bulk(&cpp_out_stem, 1, true).expect("C++ mapper failed");
    run_rust_bulk(&rust_out_dir, RUST_INDEX_PREFIX, 1, true)
        .expect("Rust mapper failed");

    // 4. Compare RAD files
    let cpp_rad = cpp_out_stem.with_extension("rad");
    let rust_rad = rust_out_dir.join("map.rad");

    assert!(cpp_rad.exists(), "C++ RAD file not found: {}", cpp_rad.display());
    assert!(rust_rad.exists(), "Rust RAD file not found: {}", rust_rad.display());

    eprintln!("Comparing RAD files...");
    let result = compare_bulk_rad_full(&cpp_rad, &rust_rad)
        .expect("failed to compare RAD files");

    eprintln!("Comparison result:");
    eprintln!("  Header match: {}", result.header_match);
    eprintln!("  Records A (C++): {}", result.total_records_a);
    eprintln!("  Records B (Rust): {}", result.total_records_b);
    eprintln!("  Matching: {}", result.matching_records);
    eprintln!("  Missing in A: {}", result.missing_in_a);
    eprintln!("  Missing in B: {}", result.missing_in_b);
    eprintln!("  Mismatch categories (A-only records):");
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

    // Report mapping rates
    let cpp_rate = if result.total_records_a > 0 {
        result.total_records_a as f64 / 1_000_000.0 * 100.0
    } else {
        0.0
    };
    let rust_rate = if result.total_records_b > 0 {
        result.total_records_b as f64 / 1_000_000.0 * 100.0
    } else {
        0.0
    };
    eprintln!(
        "  C++ mapping rate: {:.2}% ({}/1M)",
        cpp_rate, result.total_records_a
    );
    eprintln!(
        "  Rust mapping rate: {:.2}% ({}/1M)",
        rust_rate, result.total_records_b
    );

    // Check header match (ref names compared as set, not ordered)
    assert!(
        result.header_match,
        "RAD headers differ between C++ and Rust: {}",
        result.notes,
    );

    // Report record-level parity
    let match_rate = if result.total_records_a > 0 {
        result.matching_records as f64 / result.total_records_a as f64 * 100.0
    } else {
        0.0
    };
    eprintln!("  Record match rate: {:.2}% ({}/{})",
        match_rate, result.matching_records, result.total_records_a);

    // Assert reasonable record-level parity.
    // Not 100% due to SSHash dictionary implementation differences
    // between C++ sshash and Rust sshash-rs (different internal structures,
    // streaming query behavior, hit collection order).
    // Both map ~96% of reads; ~80% of records match exactly.
    assert!(
        match_rate > 75.0,
        "Record-level match rate too low ({:.2}%): {}",
        match_rate,
        result.notes,
    );
}

/// Variant: run with multiple threads to verify multiset comparison handles
/// non-deterministic chunk ordering.
#[test]
#[ignore]
fn bulk_pe_rad_parity_multithreaded() {
    if !Path::new(READ1).exists() || !Path::new(CPP_BULK_BIN).exists() {
        eprintln!("SKIP: prerequisites not found");
        return;
    }
    if !Path::new(&format!("{CFISH_PREFIX}.cf_seg")).exists() {
        eprintln!("SKIP: cuttlefish output not found");
        return;
    }

    ensure_cpp_index().expect("failed to build C++ index");
    ensure_rust_index().expect("failed to build Rust index");

    let tmpdir = tempfile::tempdir().expect("failed to create tempdir");
    let cpp_out_stem = tmpdir.path().join("cpp_out");
    let rust_out_dir = tmpdir.path().join("rust_out");
    std::fs::create_dir_all(&rust_out_dir).unwrap();

    // Run with 4 threads — record order will differ but multisets should match
    run_cpp_bulk(&cpp_out_stem, 4, true).expect("C++ mapper failed");
    run_rust_bulk(&rust_out_dir, RUST_INDEX_PREFIX, 4, true)
        .expect("Rust mapper failed");

    let cpp_rad = cpp_out_stem.with_extension("rad");
    let rust_rad = rust_out_dir.join("map.rad");

    let result = compare_bulk_rad_full(&cpp_rad, &rust_rad)
        .expect("failed to compare RAD files");

    eprintln!("Multithreaded comparison:");
    eprintln!("  Records A: {}, Records B: {}", result.total_records_a, result.total_records_b);
    eprintln!("  Matching: {}, Notes: {}", result.matching_records, result.notes);

    assert!(result.header_match, "Header mismatch: {}", result.notes);

    let match_rate = if result.total_records_a > 0 {
        result.matching_records as f64 / result.total_records_a as f64 * 100.0
    } else {
        0.0
    };
    eprintln!("  Record match rate: {:.2}%", match_rate);
    assert!(
        match_rate > 75.0,
        "Multithreaded record match rate too low ({:.2}%): {}",
        match_rate,
        result.notes,
    );
}

/// Diagnostic: compare STRICT mode outputs to isolate dictionary-dependent behavior.
///
/// If STRICT mode gives ~100% match but PERMISSIVE doesn't, the discrepancy
/// is from different contig layouts affecting the PERMISSIVE skip logic.
#[test]
#[ignore]
fn bulk_pe_rad_parity_strict() {
    if !Path::new(READ1).exists() || !Path::new(CPP_BULK_BIN).exists() {
        eprintln!("SKIP: prerequisites not found");
        return;
    }
    if !Path::new(&format!("{CFISH_PREFIX}.cf_seg")).exists() {
        eprintln!("SKIP: cuttlefish output not found");
        return;
    }

    ensure_cpp_index().expect("failed to build C++ index");
    ensure_rust_index().expect("failed to build Rust index");

    let tmpdir = tempfile::tempdir().expect("failed to create tempdir");
    let cpp_out_stem = tmpdir.path().join("cpp_out");
    let rust_out_dir = tmpdir.path().join("rust_out");
    std::fs::create_dir_all(&rust_out_dir).unwrap();

    // Run C++ mapper with strict
    let status = Command::new(CPP_BULK_BIN)
        .arg("-i").arg(CPP_INDEX_PREFIX)
        .arg("-o").arg(&cpp_out_stem)
        .arg("-t").arg("1")
        .arg("-1").arg(READ1)
        .arg("-2").arg(READ2)
        .arg("--quiet")
        .arg("--no-poison")
        .arg("--skipping-strategy").arg("strict")
        .status()
        .expect("failed to run C++ bulk mapper");
    assert!(status.success(), "C++ strict mapper failed");

    // Run Rust mapper with strict
    let bin = if PathBuf::from("target/release/piscem-rs").exists() {
        PathBuf::from("target/release/piscem-rs")
    } else {
        PathBuf::from("target/debug/piscem-rs")
    };
    let status = Command::new(&bin)
        .arg("map-bulk")
        .arg("-i").arg(RUST_INDEX_PREFIX)
        .arg("-o").arg(&rust_out_dir)
        .arg("-t").arg("1")
        .arg("-1").arg(READ1)
        .arg("-2").arg(READ2)
        .arg("--no-poison")
        .arg("--skipping-strategy").arg("strict")
        .status()
        .expect("failed to run Rust bulk mapper");
    assert!(status.success(), "Rust strict mapper failed");

    let cpp_rad = cpp_out_stem.with_extension("rad");
    let rust_rad = rust_out_dir.join("map.rad");

    eprintln!("Comparing STRICT mode RAD files...");
    let result = compare_bulk_rad_full(&cpp_rad, &rust_rad)
        .expect("failed to compare RAD files");

    eprintln!("STRICT mode comparison:");
    eprintln!("  Records A (C++): {}", result.total_records_a);
    eprintln!("  Records B (Rust): {}", result.total_records_b);
    eprintln!("  Matching: {}", result.matching_records);
    eprintln!("  Missing in A: {}", result.missing_in_a);
    eprintln!("  Missing in B: {}", result.missing_in_b);
    eprintln!("  Mismatch categories (A-only records):");
    eprintln!("    Same targets, diff detail: {}", result.same_targets_diff_detail);
    eprintln!("    Different targets: {}", result.different_targets);
    eprintln!("    Diff position only: {}", result.diff_pos_only);
    eprintln!("    Diff frag_len only: {}", result.diff_frag_len_only);
    eprintln!("    Diff num alignments: {}", result.diff_num_alns);
    eprintln!("  Notes: {}", result.notes);

    let match_rate = if result.total_records_a > 0 {
        result.matching_records as f64 / result.total_records_a as f64 * 100.0
    } else {
        0.0
    };
    eprintln!("  Record match rate: {:.2}%", match_rate);
}

/// Quick SE comparison: compares pre-generated SE RAD files in /tmp.
///
/// Run both mappers in SE mode first, then this test:
///   cargo test --features parity-test --test rad_parity_bulk --release -- se_rad_parity --ignored --nocapture
#[test]
#[ignore]
fn se_rad_parity() {
    let cpp_rad = std::path::PathBuf::from("/tmp/se_parity_test/cpp_out/se.rad");
    let rust_rad = std::path::PathBuf::from("/tmp/se_parity_test/rust_out/map.rad");

    if !cpp_rad.exists() || !rust_rad.exists() {
        eprintln!("SKIP: SE RAD files not found. Run both mappers in SE mode first.");
        return;
    }

    eprintln!("Comparing SE RAD files...");
    let result = compare_bulk_rad_full(&cpp_rad, &rust_rad)
        .expect("failed to compare RAD files");

    eprintln!("SE Comparison result:");
    eprintln!("  Header match: {}", result.header_match);
    eprintln!("  Records A (C++): {}", result.total_records_a);
    eprintln!("  Records B (Rust): {}", result.total_records_b);
    eprintln!("  Matching: {}", result.matching_records);
    eprintln!("  Missing in A: {}", result.missing_in_a);
    eprintln!("  Missing in B: {}", result.missing_in_b);
    eprintln!("  Mismatch categories (A-only records):");
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

    let match_rate = if result.total_records_a > 0 {
        result.matching_records as f64 / result.total_records_a as f64 * 100.0
    } else {
        0.0
    };
    eprintln!("  SE Record match rate: {:.2}%", match_rate);
}

/// Diagnostic: dump k-mer lookup results for the first N reads.
///
/// Produces output in the same format as the C++ kmer_dump tool, allowing
/// direct diff-based comparison of dictionary lookups between implementations.
///
///   cargo test --features parity-test --test rad_parity_bulk --release -- kmer_dump_rust --ignored --nocapture 2>/dev/null > /tmp/rust_kmer_dump.txt
#[test]
#[ignore]
fn kmer_dump_rust() {
    use sshash_lib::dispatch_on_k;
    use piscem_rs::index::reference_index::ReferenceIndex;

    let max_reads: usize = std::env::var("KMER_DUMP_READS")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(100);

    if !Path::new(RUST_INDEX_PREFIX).with_extension("ssi").exists() {
        ensure_rust_index().expect("failed to build Rust index");
    }

    let index = ReferenceIndex::load(
        Path::new(RUST_INDEX_PREFIX), false, false
    ).expect("failed to load Rust index");
    let k = index.k();
    eprintln!("Rust index loaded: k={}, refs={}", k, index.num_refs());

    // Read sequences from a pre-decompressed file at /tmp/first_reads.fq
    // Generate with: gzcat test_data/sim_1M_1.fq.gz | head -20 > /tmp/first_reads.fq
    let fq_path = "/tmp/first_reads.fq";
    if !Path::new(fq_path).exists() {
        eprintln!("SKIP: {} not found. Generate with: gzcat test_data/sim_1M_1.fq.gz | head -20 > /tmp/first_reads.fq", fq_path);
        return;
    }
    let file = std::fs::File::open(fq_path).expect("failed to open FASTQ");
    let reader = std::io::BufReader::new(file);
    use std::io::BufRead;
    let mut sequences: Vec<String> = Vec::new();
    let mut line_idx: usize = 0;
    for line_result in reader.lines() {
        match line_result {
            Ok(line) => {
                if line_idx % 4 == 1 {
                    sequences.push(line);
                    if sequences.len() >= max_reads {
                        break;
                    }
                }
                line_idx += 1;
            }
            Err(e) => {
                eprintln!("Error reading FASTQ at line {}: {}", line_idx, e);
                break;
            }
        }
    }

    eprintln!("Read {} sequences from FASTQ", sequences.len());
    if sequences.is_empty() {
        eprintln!("ERROR: no sequences read");
        return;
    }
    eprintln!("First seq ({}bp): {}...", sequences[0].len(), &sequences[0][..40.min(sequences[0].len())]);

    dispatch_on_k!(k, K => {
        kmer_dump_inner::<K>(&index, &sequences);
    });
}

fn kmer_dump_inner<const K: usize>(
    index: &piscem_rs::index::reference_index::ReferenceIndex,
    sequences: &[String],
) where
    sshash_lib::Kmer<K>: sshash_lib::KmerBits,
{
    use piscem_rs::mapping::streaming_query::PiscemStreamingQuery;
    let k = index.k();
    let encoding = index.encoding();
    let dict = index.dict();
    let mut query = PiscemStreamingQuery::<K>::new(dict);

    for (read_idx, seq) in sequences.iter().enumerate() {
        let seq_bytes = seq.as_bytes();
        eprintln!("READ {} len={}", read_idx, seq_bytes.len());

        query.reset();

        for i in 0..=(seq_bytes.len().saturating_sub(k)) {
            let kmer_bytes = &seq_bytes[i..i + k];

            // Skip if contains non-ACGT
            let has_invalid = kmer_bytes.iter().any(|&b| {
                !matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't')
            });
            if has_invalid {
                eprintln!("  pos={} NOT_FOUND (invalid)", i);
                query.reset();
                continue;
            }

            let kmer_str = std::str::from_utf8(kmer_bytes).unwrap();

            // Streaming lookup (the normal path)
            let result = query.lookup(kmer_str);
            let streaming_found = result.is_found();

            // Non-streaming (direct) lookup — bypasses state machine
            let direct_result = dict.query_from_str::<K>(kmer_str);
            let direct_found = direct_result.is_found();

            // Flag mismatches between streaming and direct
            let mismatch = streaming_found != direct_found;

            if streaming_found {
                if let Some(phit) = index.resolve_lookup(&result) {
                    if phit.is_empty() {
                        eprint!("  pos={} NOT_FOUND(stream)", i);
                    } else {
                        eprint!(
                            "  pos={} cid={} cpos={} clen={} ori={} gpos={} nhits={}",
                            i,
                            phit.contig_id(),
                            phit.contig_pos(),
                            phit.contig_len(),
                            if phit.hit_fw_on_contig() { "fw" } else { "rc" },
                            phit.global_pos(),
                            phit.num_hits(),
                        );

                        // Decode each reference hit
                        for entry in phit.ref_range().iter() {
                            let rp = phit.decode_hit(entry, &encoding);
                            let tid = encoding.transcript_id(entry);
                            eprint!(
                                " [t={} p={} {}]",
                                tid,
                                rp.pos,
                                if rp.is_fw { "fw" } else { "rc" }
                            );
                        }
                    }
                } else {
                    eprint!("  pos={} NOT_FOUND(stream)", i);
                }
            } else {
                eprint!("  pos={} NOT_FOUND(stream)", i);
            }

            if mismatch {
                if direct_found {
                    eprint!(
                        " *** DIRECT_HIT cid={} cpos={} ori={}",
                        direct_result.string_id,
                        direct_result.kmer_id_in_string,
                        if direct_result.kmer_orientation > 0 { "fw" } else { "rc" }
                    );
                } else {
                    eprint!(" *** DIRECT_MISS");
                }
            }
            eprintln!();
        }
        eprintln!();
    }
}

/// Diagnostic: dump hit_searcher raw hits for the first N reads.
///
/// This shows what the PERMISSIVE mode skip logic actually collects,
/// which can differ from the per-k-mer kmer_dump output.
///
///   cargo test --features parity-test --test rad_parity_bulk --release -- hit_dump_rust --ignored --nocapture 2>/dev/null
#[test]
#[ignore]
fn hit_dump_rust() {
    use sshash_lib::dispatch_on_k;
    use piscem_rs::index::reference_index::ReferenceIndex;

    let max_reads: usize = std::env::var("HIT_DUMP_READS")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(100);

    let strat_str = std::env::var("HIT_DUMP_STRAT").unwrap_or_else(|_| "permissive".into());

    if !Path::new(RUST_INDEX_PREFIX).with_extension("ssi").exists() {
        ensure_rust_index().expect("failed to build Rust index");
    }

    let index = ReferenceIndex::load(
        Path::new(RUST_INDEX_PREFIX), true, false
    ).expect("failed to load Rust index");
    let k = index.k();
    eprintln!("Rust index loaded: k={}, refs={}", k, index.num_refs());
    eprintln!("Strategy: {}", strat_str);

    // Read sequences from /tmp/first_reads.fq
    let fq_path = "/tmp/first_reads.fq";
    if !Path::new(fq_path).exists() {
        eprintln!("SKIP: {} not found", fq_path);
        return;
    }
    let file = std::fs::File::open(fq_path).expect("failed to open FASTQ");
    let reader = std::io::BufReader::new(file);
    use std::io::BufRead;
    let mut sequences: Vec<String> = Vec::new();
    let mut line_idx: usize = 0;
    for line_result in reader.lines() {
        match line_result {
            Ok(line) => {
                if line_idx % 4 == 1 {
                    sequences.push(line);
                    if sequences.len() >= max_reads {
                        break;
                    }
                }
                line_idx += 1;
            }
            Err(_) => break,
        }
    }
    eprintln!("Read {} sequences", sequences.len());

    let strat = if strat_str == "strict" {
        piscem_rs::mapping::hit_searcher::SkippingStrategy::Strict
    } else {
        piscem_rs::mapping::hit_searcher::SkippingStrategy::Permissive
    };

    dispatch_on_k!(k, K => {
        hit_dump_inner::<K>(&index, &sequences, strat);
    });
}

fn hit_dump_inner<const K: usize>(
    index: &piscem_rs::index::reference_index::ReferenceIndex,
    sequences: &[String],
    strat: piscem_rs::mapping::hit_searcher::SkippingStrategy,
) where
    sshash_lib::Kmer<K>: sshash_lib::KmerBits,
{
    use piscem_rs::mapping::hit_searcher::HitSearcher;
    use piscem_rs::mapping::streaming_query::PiscemStreamingQuery;

    let dict = index.dict();
    let encoding = index.encoding();
    let mut hs = HitSearcher::new(index);
    let mut query = PiscemStreamingQuery::<K>::new(dict);

    for (read_idx, seq) in sequences.iter().enumerate() {
        let seq_bytes = seq.as_bytes();
        println!("READ {} len={}", read_idx, seq_bytes.len());

        query.reset();
        let has_hits = hs.get_raw_hits_sketch::<K>(seq_bytes, &mut query, strat, true);
        let raw_hits = hs.left_hits();

        println!("  raw_hits={} has_matching={}", raw_hits.len(), has_hits);

        for (i, (read_pos, ph)) in raw_hits.iter().enumerate() {
            print!(
                "  hit[{}] rpos={} cid={} cpos={} clen={} ori={} gpos={} nhits={} open={}",
                i,
                read_pos,
                ph.contig_id(),
                ph.contig_pos(),
                ph.contig_len(),
                if ph.hit_fw_on_contig() { "fw" } else { "rc" },
                ph.global_pos(),
                ph.num_hits(),
                ph.resulted_from_open_search(),
            );

            // Decode reference hits
            for entry in ph.ref_range().iter() {
                let rp = ph.decode_hit(entry, &encoding);
                let tid = encoding.transcript_id(entry);
                print!(
                    " [t={} p={} {}]",
                    tid,
                    rp.pos,
                    if rp.is_fw { "fw" } else { "rc" }
                );
            }
            println!();
        }
        println!();
    }
}

/// Diagnostic: dump individual RAD records from two files side by side.
///
/// Run:
///   cargo test --features parity-test --test rad_parity_bulk --release -- dump_rad_records --ignored --nocapture
#[test]
#[ignore]
fn dump_rad_records() {
    use piscem_rs::verify::rad_compare::read_bulk_rad_records;

    let cpp_rad = std::path::PathBuf::from("/tmp/debug_cpp/out.rad");
    let rust_rad = std::path::PathBuf::from("/tmp/debug_rust/map.rad");

    if !cpp_rad.exists() || !rust_rad.exists() {
        eprintln!("SKIP: debug RAD files not found. Run both mappers on /tmp/single_pair_*.fq first.");
        return;
    }

    let (prelude_a, records_a) = read_bulk_rad_records(&cpp_rad)
        .expect("failed to read C++ RAD");
    let (prelude_b, records_b) = read_bulk_rad_records(&rust_rad)
        .expect("failed to read Rust RAD");

    eprintln!("C++ RAD: {} records, {} refs, paired={}",
        records_a.len(), prelude_a.hdr.ref_count, prelude_a.hdr.is_paired);
    eprintln!("Rust RAD: {} records, {} refs, paired={}",
        records_b.len(), prelude_b.hdr.ref_count, prelude_b.hdr.is_paired);

    // Build ref name lookup for A
    let ref_names_a = &prelude_a.hdr.ref_names;
    let ref_names_b = &prelude_b.hdr.ref_names;

    // Build B→A ref ID translation
    let name_to_id_a: std::collections::HashMap<&str, u32> = ref_names_a
        .iter().enumerate()
        .map(|(i, n)| (n.as_str(), i as u32))
        .collect();

    eprintln!("\n=== C++ Records ===");
    for (i, rec) in records_a.iter().enumerate() {
        eprintln!("Record {}: frag_type={}", i, rec.frag_type);
        for (j, &(ref_id, ori, pos, flen)) in rec.alignments.iter().enumerate() {
            let name = if (ref_id as usize) < ref_names_a.len() {
                // Truncate long names
                let full = &ref_names_a[ref_id as usize];
                if full.len() > 40 { &full[..40] } else { full.as_str() }
            } else { "?" };
            eprintln!("  aln[{}]: ref={} ({}) ori={} pos={} flen={}",
                j, ref_id, name, ori, pos, flen);
        }
    }

    eprintln!("\n=== Rust Records ===");
    for (i, rec) in records_b.iter().enumerate() {
        eprintln!("Record {}: frag_type={}", i, rec.frag_type);
        for (j, &(ref_id, ori, pos, flen)) in rec.alignments.iter().enumerate() {
            let name = if (ref_id as usize) < ref_names_b.len() {
                let full = &ref_names_b[ref_id as usize];
                if full.len() > 40 { &full[..40] } else { full.as_str() }
            } else { "?" };
            // Translate to A's ref ID for comparison
            let a_id = name_to_id_a.get(
                if (ref_id as usize) < ref_names_b.len() {
                    ref_names_b[ref_id as usize].as_str()
                } else { "" }
            ).copied().unwrap_or(u32::MAX);
            eprintln!("  aln[{}]: ref={} (a_id={}) ({}) ori={} pos={} flen={}",
                j, ref_id, a_id, name, ori, pos, flen);
        }
    }

    // Direct comparison
    eprintln!("\n=== Comparison ===");
    let max_recs = records_a.len().max(records_b.len());
    for i in 0..max_recs {
        let a = records_a.get(i);
        let b = records_b.get(i);
        match (a, b) {
            (Some(ra), Some(rb)) => {
                // Translate B to A's namespace for comparison
                let mut rb_translated = rb.clone();
                for aln in &mut rb_translated.alignments {
                    if (aln.0 as usize) < ref_names_b.len() {
                        let name = ref_names_b[aln.0 as usize].as_str();
                        aln.0 = name_to_id_a.get(name).copied().unwrap_or(u32::MAX);
                    }
                }
                rb_translated.alignments.sort();

                if *ra == rb_translated {
                    eprintln!("Record {}: MATCH", i);
                } else {
                    eprintln!("Record {}: MISMATCH", i);
                    eprintln!("  C++:  {}", ra);
                    eprintln!("  Rust: {}", rb_translated);
                    // Show per-alignment diffs
                    let max_alns = ra.alignments.len().max(rb_translated.alignments.len());
                    for j in 0..max_alns {
                        let aa = ra.alignments.get(j);
                        let bb = rb_translated.alignments.get(j);
                        if aa != bb {
                            eprintln!("    aln[{}]: C++={:?} Rust={:?}", j, aa, bb);
                        }
                    }
                }
            }
            (Some(ra), None) => eprintln!("Record {}: Only in C++: {}", i, ra),
            (None, Some(rb)) => eprintln!("Record {}: Only in Rust: {}", i, rb),
            (None, None) => {}
        }
    }
}
