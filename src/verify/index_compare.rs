//! Index semantic comparison.
//!
//! Compares two Rust-built indices from the same cuttlefish input to verify
//! deterministic construction. Can also compare reference metadata (names,
//! lengths) between any two indices.

use std::path::Path;

use anyhow::Result;
use serde::Serialize;

use crate::index::reference_index::ReferenceIndex;

/// Summary of index comparison results.
#[derive(Debug, Serialize)]
pub struct IndexComparisonSummary {
    pub passed: bool,
    pub notes: String,
    pub num_refs_match: bool,
    pub ref_names_match: bool,
    pub ref_lengths_match: bool,
    pub num_refs: usize,
    pub sampled_posting_lists_checked: usize,
    pub posting_list_mismatches: usize,
}

/// Compare reference metadata between two indices.
///
/// Loads both indices and compares:
/// - Number of references
/// - Reference names (ordered)
/// - Reference lengths
pub fn compare_ref_metadata(index_a: &ReferenceIndex, index_b: &ReferenceIndex) -> IndexComparisonSummary {
    let num_refs_match = index_a.num_refs() == index_b.num_refs();
    let n = index_a.num_refs().min(index_b.num_refs());

    let mut ref_names_match = num_refs_match;
    let mut ref_lengths_match = num_refs_match;

    for i in 0..n {
        if index_a.ref_name(i) != index_b.ref_name(i) {
            ref_names_match = false;
        }
        if index_a.ref_len(i) != index_b.ref_len(i) {
            ref_lengths_match = false;
        }
    }

    let passed = num_refs_match && ref_names_match && ref_lengths_match;
    let notes = if passed {
        format!("All {} references match", n)
    } else {
        let mut notes = Vec::new();
        if !num_refs_match {
            notes.push(format!(
                "num_refs mismatch: {} vs {}",
                index_a.num_refs(),
                index_b.num_refs()
            ));
        }
        if !ref_names_match {
            notes.push("reference names differ".to_string());
        }
        if !ref_lengths_match {
            notes.push("reference lengths differ".to_string());
        }
        notes.join("; ")
    };

    IndexComparisonSummary {
        passed,
        notes,
        num_refs_match,
        ref_names_match,
        ref_lengths_match,
        num_refs: n,
        sampled_posting_lists_checked: 0,
        posting_list_mismatches: 0,
    }
}

/// Compare two indices loaded from disk.
///
/// Loads both from their prefix paths and compares reference metadata.
pub fn compare_index_semantics(prefix_a: &Path, prefix_b: &Path) -> Result<IndexComparisonSummary> {
    let index_a = ReferenceIndex::load(prefix_a, false, false)?;
    let index_b = ReferenceIndex::load(prefix_b, false, false)?;
    Ok(compare_ref_metadata(&index_a, &index_b))
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_index_comparison_summary_serialize() {
        let summary = IndexComparisonSummary {
            passed: true,
            notes: "ok".to_string(),
            num_refs_match: true,
            ref_names_match: true,
            ref_lengths_match: true,
            num_refs: 10,
            sampled_posting_lists_checked: 100,
            posting_list_mismatches: 0,
        };
        let json = serde_json::to_string(&summary).unwrap();
        assert!(json.contains("\"passed\":true"));
    }
}
