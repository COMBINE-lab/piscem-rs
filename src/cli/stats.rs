use anyhow::{Context, Result};
use clap::Args;
use serde_json::{json, Value};
use std::path::Path;

use crate::index::reference_index::ReferenceIndex;

#[derive(Args, Debug)]
#[command(about = "Print index statistics as pretty-printed JSON")]
pub struct StatsArgs {
    /// Path prefix for the index files (e.g. /path/to/index)
    #[arg(short, long)]
    pub index: String,
}

/// File extensions that make up the index on disk.
const INDEX_EXTENSIONS: &[(&str, &str)] = &[
    ("ssi", "dictionary"),
    ("ssi.mphf", "mphf"),
    ("ctab", "contig_table"),
    ("refinfo", "ref_info"),
    ("ectab", "eq_class_table"),
    ("poison", "poison_table"),
    ("poison.json", "poison_stats"),
    ("sigs.json", "signatures"),
];

fn file_size_bytes(prefix: &Path, ext: &str) -> Option<u64> {
    let mut p = prefix.to_path_buf();
    p.set_extension(ext);
    std::fs::metadata(&p).ok().map(|m| m.len())
}

pub fn run(args: StatsArgs) -> Result<()> {
    let prefix = Path::new(&args.index);

    // Load the full index (including EC table and poison table if present).
    let index = ReferenceIndex::load(prefix, true, true)
        .with_context(|| format!("failed to load index at {}", prefix.display()))?;

    let dict = index.dict();
    let ctab = index.contig_table();

    // Reference length stats.
    let num_refs = index.num_refs();
    let total_ref_len: u64 = (0..num_refs).map(|i| index.ref_len(i)).sum();
    let max_ref_len = index.ref_info().max_ref_len();
    let min_ref_len = (0..num_refs).map(|i| index.ref_len(i)).min().unwrap_or(0);
    let mean_ref_len = if num_refs > 0 {
        total_ref_len as f64 / num_refs as f64
    } else {
        0.0
    };

    // Contig table scan.
    let num_contigs = ctab.num_contigs();
    let num_entries = ctab.num_entries();
    let mut max_entries_per_contig: usize = 0;
    let mut contigs_with_zero_entries: usize = 0;
    for c in 0..num_contigs {
        let count = ctab.contig_entries(c as u64).len();
        if count > max_entries_per_contig {
            max_entries_per_contig = count;
        }
        if count == 0 {
            contigs_with_zero_entries += 1;
        }
    }

    // On-disk file sizes.
    let mut file_sizes = serde_json::Map::new();
    let mut total_disk_bytes: u64 = 0;
    for &(ext, label) in INDEX_EXTENSIONS {
        if let Some(size) = file_size_bytes(prefix, ext) {
            file_sizes.insert(label.to_string(), json!(size));
            total_disk_bytes += size;
        }
    }
    file_sizes.insert("total".to_string(), json!(total_disk_bytes));

    // In-memory sizes.
    let mut mem_sizes = serde_json::Map::new();
    mem_sizes.insert("contig_table".to_string(), json!(ctab.size_bytes()));
    if let Some(ec) = index.ec_table() {
        mem_sizes.insert("eq_class_table".to_string(), json!(ec.size_bytes()));
    }

    // Build the output object.
    let encoding = index.encoding();
    let entry_width_bits = ctab.ref_len_bits() + ctab.num_ref_bits() + 1;

    let mut root = serde_json::Map::new();

    root.insert("k".to_string(), json!(index.k()));
    root.insert("m".to_string(), json!(index.m()));

    root.insert("references".to_string(), json!({
        "num_refs": num_refs,
        "total_length": total_ref_len,
        "min_length": min_ref_len,
        "max_length": max_ref_len,
        "mean_length": (mean_ref_len * 100.0).round() / 100.0,
    }));

    root.insert("dictionary".to_string(), json!({
        "num_unitigs": dict.num_strings(),
        "num_minimizers": dict.num_minimizers(),
        "num_bits": dict.num_bits(),
        "canonical": dict.canonical(),
    }));

    root.insert("contig_table".to_string(), json!({
        "num_contigs": num_contigs,
        "num_entries": num_entries,
        "entry_width_bits": entry_width_bits,
        "ref_id_bits": ctab.num_ref_bits(),
        "position_bits": ctab.ref_len_bits(),
        "ref_shift": encoding.ref_shift,
        "pos_mask": format!("0x{:x}", encoding.pos_mask),
        "max_entries_per_contig": max_entries_per_contig,
        "contigs_with_zero_entries": contigs_with_zero_entries,
    }));

    if let Some(ec) = index.ec_table() {
        root.insert("eq_class_table".to_string(), json!({
            "num_tiles": ec.num_tiles(),
            "num_ecs": ec.num_ecs(),
            "num_label_entries": ec.num_label_entries(),
        }));
    }

    if let Some(pt) = index.poison_table() {
        root.insert("poison_table".to_string(), json!({
            "num_poison_kmers": pt.num_poison_kmers(),
            "num_poison_occurrences": pt.num_poison_occs(),
        }));
    }

    root.insert("size_on_disk_bytes".to_string(), Value::Object(file_sizes));
    root.insert("size_in_memory_bytes".to_string(), Value::Object(mem_sizes));

    let output = serde_json::to_string_pretty(&Value::Object(root))?;
    println!("{output}");

    Ok(())
}
