use std::fs;

use anyhow::Result;
use serde::Serialize;

use crate::verify::index_compare::{compare_index_semantics, IndexComparisonSummary};
use crate::verify::rad_compare::{compare_rad_semantics, RadComparisonSummary};

#[derive(Debug)]
pub struct Phase0ParityConfig {
    pub dataset: String,
    pub cpp_bin_dir: Option<String>,
    pub output_report: Option<String>,
}

#[derive(Debug, Serialize)]
pub struct Phase0ParityReport {
    pub dataset: String,
    pub cpp_bin_dir: Option<String>,
    pub index_summary: IndexComparisonSummary,
    pub rad_summary: RadComparisonSummary,
    pub passed: bool,
}

pub fn run_phase0_parity(cfg: Phase0ParityConfig) -> Result<()> {
    let index_summary = compare_index_semantics(&cfg.dataset)?;

    let cpp_rad = format!("{}/cpp.rad", cfg.dataset);
    let rust_rad = format!("{}/rust.rad", cfg.dataset);
    let rad_summary = compare_rad_semantics(&cpp_rad, &rust_rad)?;

    let passed = index_summary.passed && rad_summary.passed;
    let report = Phase0ParityReport {
        dataset: cfg.dataset,
        cpp_bin_dir: cfg.cpp_bin_dir,
        index_summary,
        rad_summary,
        passed,
    };

    let report_path = cfg
        .output_report
        .unwrap_or_else(|| "phase0_parity_report.json".to_string());

    fs::write(&report_path, serde_json::to_vec_pretty(&report)?)?;

    if report.passed {
        tracing::info!(report_path, "Phase 0 parity check passed");
        Ok(())
    } else {
        anyhow::bail!(
            "Phase 0 parity check failed. See report: {}",
            report_path
        )
    }
}
