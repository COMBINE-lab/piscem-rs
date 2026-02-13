use anyhow::Result;
use serde::Serialize;

#[derive(Debug, Serialize)]
pub struct RadComparisonSummary {
    pub passed: bool,
    pub notes: String,
}

pub fn compare_rad_semantics(_cpp_rad: &str, _rust_rad: &str) -> Result<RadComparisonSummary> {
    // Phase 0 stub: RAD comparison will be implemented when libradicl
    // integration is wired up.
    Ok(RadComparisonSummary {
        passed: false,
        notes: "Phase 0 stub: libradicl RAD decoder not yet wired".to_string(),
    })
}
