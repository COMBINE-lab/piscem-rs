//! map_info.json writer — summary statistics for a mapping run.

use std::path::Path;

use anyhow::Result;

use crate::index::reference_index::RefSigInfo;

/// All parameters needed to write a map_info.json file.
pub struct MapInfoParams<'a> {
    pub path: &'a Path,
    /// Mapping mode: "bulk", "sc-rna", or "sc-atac".
    pub mode: &'a str,
    pub num_reads: u64,
    pub num_mapped: u64,
    pub num_poisoned: u64,
    pub elapsed_secs: f64,
    pub sig_info: Option<&'a RefSigInfo>,
    pub piscem_rs_version: &'a str,
    pub num_threads: usize,
    pub index_path: &'a Path,
    pub k: usize,
    pub m: usize,
    pub num_refs: usize,
    /// Skipping strategy: "permissive", "strict", or "every-kmer".
    pub skipping_strategy: &'a str,
}

/// Write a map_info.json file with mapping statistics and run metadata.
pub fn write_map_info(params: &MapInfoParams) -> Result<()> {
    let percent_mapped = if params.num_reads > 0 {
        (params.num_mapped as f64 / params.num_reads as f64) * 100.0
    } else {
        0.0
    };

    let cmdline: Vec<String> = std::env::args().collect();

    let mut info = serde_json::json!({
        "mode": params.mode,
        "piscem_rs_version": params.piscem_rs_version,
        "index_path": params.index_path.display().to_string(),
        "k": params.k,
        "m": params.m,
        "num_refs": params.num_refs,
        "num_threads": params.num_threads,
        "skipping_strategy": params.skipping_strategy,
        "num_reads": params.num_reads,
        "num_mapped": params.num_mapped,
        "num_poisoned": params.num_poisoned,
        "percent_mapped": format!("{:.2}", percent_mapped),
        "runtime_seconds": format!("{:.2}", params.elapsed_secs),
        "cmdline": cmdline.join(" "),
    });

    if let Some(sig) = params.sig_info {
        info["signatures"] = serde_json::json!({
            "sha256_names": sig.sha256_names,
            "sha256_seqs":  sig.sha256_seqs,
            "sha512_names": sig.sha512_names,
            "sha512_seqs":  sig.sha512_seqs,
        });
    }

    let file = std::fs::File::create(params.path)?;
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
        let index_path = std::path::PathBuf::from("/tmp/test_index");

        write_map_info(&MapInfoParams {
            path: &path,
            mode: "bulk",
            num_reads: 1000,
            num_mapped: 800,
            num_poisoned: 10,
            elapsed_secs: 42.5,
            sig_info: None,
            piscem_rs_version: "0.1.0",
            num_threads: 8,
            index_path: &index_path,
            k: 31,
            m: 19,
            num_refs: 500,
            skipping_strategy: "permissive",
        })
        .unwrap();

        let content = std::fs::read_to_string(&path).unwrap();
        let val: serde_json::Value = serde_json::from_str(&content).unwrap();
        assert_eq!(val["num_reads"], 1000);
        assert_eq!(val["num_mapped"], 800);
        assert_eq!(val["num_poisoned"], 10);
        assert_eq!(val["percent_mapped"], "80.00");
        assert_eq!(val["mode"], "bulk");
        assert_eq!(val["k"], 31);
        assert_eq!(val["m"], 19);
        assert_eq!(val["num_refs"], 500);
        assert_eq!(val["num_threads"], 8);
        assert_eq!(val["skipping_strategy"], "permissive");
        assert_eq!(val["piscem_rs_version"], "0.1.0");
    }
}
