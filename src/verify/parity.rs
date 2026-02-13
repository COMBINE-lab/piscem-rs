//! Parity verification — end-to-end comparison of Rust outputs.
//!
//! Orchestrates index and RAD file comparisons, producing a report.

use std::fs;
use std::path::Path;

use anyhow::Result;
use serde::Serialize;

use crate::verify::index_compare::IndexComparisonSummary;
use crate::verify::rad_compare::{compare_rad_files, RadComparisonSummary};

/// Configuration for a parity check run.
#[derive(Debug)]
pub struct ParityConfig {
    /// Dataset directory or prefix.
    pub dataset: String,
    /// Optional path to C++ piscem binary directory.
    pub cpp_bin_dir: Option<String>,
    /// Path for output report JSON.
    pub output_report: Option<String>,
}

/// Full parity report.
#[derive(Debug, Serialize)]
pub struct ParityReport {
    pub dataset: String,
    pub cpp_bin_dir: Option<String>,
    pub index_summary: Option<IndexComparisonSummary>,
    pub rad_summary: Option<RadComparisonSummary>,
    pub passed: bool,
    pub notes: String,
}

/// Run a parity check.
///
/// Compares index semantics and/or RAD files depending on what's available.
pub fn run_parity(cfg: ParityConfig) -> Result<()> {
    let mut notes = Vec::new();
    let mut all_passed = true;

    // Try index comparison if two index paths are found
    let index_summary = {
        let prefix_a = Path::new(&cfg.dataset);
        if prefix_a.exists() || prefix_a.with_extension("sshash").exists() {
            // We'd need a second index to compare against.
            // For now, just note that the index is loadable.
            notes.push("index path found (comparison requires two indices)".to_string());
            None
        } else {
            notes.push("no index found at dataset path".to_string());
            None
        }
    };

    // Try RAD comparison
    let rad_summary = {
        let rad_a = Path::new(&cfg.dataset).join("map.rad");
        let rad_b_dir = cfg.cpp_bin_dir.as_ref().map(|d| Path::new(d).join("map.rad"));

        if rad_a.exists() {
            if let Some(ref rad_b) = rad_b_dir {
                if rad_b.exists() {
                    match compare_rad_files(&rad_a, rad_b) {
                        Ok(summary) => {
                            if !summary.passed {
                                all_passed = false;
                            }
                            Some(summary)
                        }
                        Err(e) => {
                            notes.push(format!("RAD comparison error: {}", e));
                            all_passed = false;
                            None
                        }
                    }
                } else {
                    notes.push(format!("C++ RAD file not found: {}", rad_b.display()));
                    None
                }
            } else {
                notes.push("Rust RAD file found, no C++ path provided".to_string());
                None
            }
        } else {
            notes.push("no Rust RAD file found".to_string());
            None
        }
    };

    let report = ParityReport {
        dataset: cfg.dataset,
        cpp_bin_dir: cfg.cpp_bin_dir,
        index_summary,
        rad_summary,
        passed: all_passed,
        notes: notes.join("; "),
    };

    let report_path = cfg
        .output_report
        .unwrap_or_else(|| "parity_report.json".to_string());

    fs::write(&report_path, serde_json::to_vec_pretty(&report)?)?;

    if report.passed {
        tracing::info!(report_path, "Parity check passed");
        Ok(())
    } else {
        anyhow::bail!(
            "Parity check failed. See report: {}",
            report_path
        )
    }
}
