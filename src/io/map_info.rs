//! map_info.json writer — summary statistics for a mapping run.

use std::path::Path;

use anyhow::Result;

use crate::index::reference_index::RefSigInfo;

/// Write a map_info.json file with mapping statistics.
pub fn write_map_info(
    path: &Path,
    num_reads: u64,
    num_mapped: u64,
    num_poisoned: u64,
    elapsed_secs: f64,
    sig_info: Option<&RefSigInfo>,
) -> Result<()> {
    let percent_mapped = if num_reads > 0 {
        (num_mapped as f64 / num_reads as f64) * 100.0
    } else {
        0.0
    };

    let cmdline: Vec<String> = std::env::args().collect();

    let mut info = serde_json::json!({
        "num_reads": num_reads,
        "num_mapped": num_mapped,
        "num_poisoned": num_poisoned,
        "percent_mapped": format!("{:.2}", percent_mapped),
        "runtime_seconds": format!("{:.2}", elapsed_secs),
        "cmdline": cmdline.join(" "),
    });

    if let Some(sig) = sig_info {
        info["signatures"] = serde_json::json!({
            "sha256_names": sig.sha256_names,
            "sha256_seqs":  sig.sha256_seqs,
            "sha512_names": sig.sha512_names,
            "sha512_seqs":  sig.sha512_seqs,
        });
    }

    let file = std::fs::File::create(path)?;
    let writer = std::io::BufWriter::new(file);
    serde_json::to_writer_pretty(writer, &info)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_map_info() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("map_info.json");

        write_map_info(&path, 1000, 800, 10, 42.5, None).unwrap();

        let content = std::fs::read_to_string(&path).unwrap();
        let val: serde_json::Value = serde_json::from_str(&content).unwrap();
        assert_eq!(val["num_reads"], 1000);
        assert_eq!(val["num_mapped"], 800);
        assert_eq!(val["num_poisoned"], 10);
        assert_eq!(val["percent_mapped"], "80.00");
    }
}
