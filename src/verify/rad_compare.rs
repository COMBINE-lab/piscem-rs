use anyhow::Result;
use serde::Serialize;

use crate::io::rad::{LibradiclAdapter, RadDecoder};

#[derive(Debug, Serialize)]
pub struct RadComparisonSummary {
    pub passed: bool,
    pub notes: String,
}

pub fn compare_rad_semantics(cpp_rad: &str, rust_rad: &str) -> Result<RadComparisonSummary> {
    let adapter = LibradiclAdapter;

    let cpp_records = adapter.decode_records(cpp_rad);
    let rust_records = adapter.decode_records(rust_rad);

    match (cpp_records, rust_records) {
        (Ok(cpp), Ok(rust)) => {
            let passed = cpp == rust;
            Ok(RadComparisonSummary {
                passed,
                notes: if passed {
                    "RAD records are semantically equivalent".to_string()
                } else {
                    format!(
                        "RAD mismatch: cpp_records={}, rust_records={}",
                        cpp.len(),
                        rust.len()
                    )
                },
            })
        }
        _ => Ok(RadComparisonSummary {
            passed: false,
            notes: "Phase 0 stub: libradicl RAD decoder not yet wired".to_string(),
        }),
    }
}
