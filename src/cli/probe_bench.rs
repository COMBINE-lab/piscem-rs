//! Diagnostic: compare the calibration probe's estimate against a whole-file
//! measurement of the same work, in the same binary.
//!
//! # Why this exists
//!
//! The probe reports a producer rate from a short prefix. Validating it against
//! a standalone harness gave 0.43-0.52 GB/s where the probe said 0.70-0.83, but
//! that harness was a separate crate built with `opt-level = 3` alone, while
//! this one builds with `lto = "thin"` and `codegen-units = 1`. A parse-heavy
//! loop can differ by half between those, so the discrepancy might have been
//! entirely a toolchain artifact rather than anything wrong with the probe.
//!
//! Measuring both here removes that confound: identical compiler settings,
//! identical `paraseq` version, one process. Whatever difference remains is
//! real, and attributable to sample size or sample position rather than to how
//! the reference was built.

use anyhow::Result;
use clap::Args;
use paraseq::fastx;
use paraseq::prelude::*;
use std::path::PathBuf;
use std::time::Instant;

#[derive(Args, Debug)]
pub struct ProbeBenchArgs {
    /// Input FASTX file (optionally compressed).
    #[arg(short = 'r', long)]
    pub reads: PathBuf,
    /// Records the probe-sized measurement should consume.
    #[arg(long, default_value = "50000")]
    pub probe_records: usize,
    /// Cap on records for the whole-file measurement (0 = no cap).
    #[arg(long, default_value = "0")]
    pub max_records: usize,
    /// Repeat the probe-sized measurement this many times, to expose variance.
    #[arg(long, default_value = "5")]
    pub reps: usize,
}

/// Producer throughput over `limit` records (0 = whole file), in GB/s of
/// sequence bytes, using exactly the loop `calibrate::probe_path` uses.
fn producer_rate(path: &PathBuf, limit: usize) -> Result<(f64, usize, u64)> {
    let (r, _) = niffler::send::from_path(path)?;
    let mut reader = fastx::Reader::new(r)?;

    // Same warm-up as the probe: one set pulled and dropped.
    let mut warm = reader.new_record_set();
    let _ = warm.fill(&mut reader);
    drop(warm);

    let mut records = 0usize;
    let mut bytes = 0u64;
    let start = Instant::now();
    loop {
        if limit > 0 && records >= limit {
            break;
        }
        let mut rs = reader.new_record_set();
        match rs.fill(&mut reader) {
            Ok(true) => {}
            _ => break,
        }
        for rec in rs.iter().filter_map(Result::ok) {
            records += 1;
            bytes += rec.seq().as_ref().len() as u64;
        }
    }
    let secs = start.elapsed().as_secs_f64();
    Ok(((bytes as f64 / 1e9) / secs.max(1e-9), records, bytes))
}

pub fn run(args: ProbeBenchArgs) -> Result<()> {
    println!("input: {}", args.reads.display());

    let mut rates = Vec::with_capacity(args.reps);
    for i in 0..args.reps {
        let (rate, recs, bytes) = producer_rate(&args.reads, args.probe_records)?;
        println!(
            "  probe-sized rep {:>2}: {:>6.3} GB/s  ({} records, {:.1} MB seq)",
            i + 1,
            rate,
            recs,
            bytes as f64 / 1e6
        );
        rates.push(rate);
    }
    rates.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let median = rates[rates.len() / 2];
    println!(
        "  probe-sized median: {median:.3} GB/s   (spread {:.3} - {:.3})",
        rates[0],
        rates[rates.len() - 1]
    );

    let (whole, recs, bytes) = producer_rate(&args.reads, args.max_records)?;
    println!(
        "  whole-file        : {whole:.3} GB/s  ({} records, {:.2} GB seq)",
        recs,
        bytes as f64 / 1e9
    );

    let ratio = median / whole;
    println!(
        "\n  probe / whole-file = {ratio:.2}x  -> probe {} supply",
        if ratio > 1.05 {
            "OVERSTATES"
        } else if ratio < 0.95 {
            "understates"
        } else {
            "matches"
        }
    );
    Ok(())
}
