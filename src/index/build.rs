//! Index build pipeline — constructs a piscem index from cuttlefish output.
//!
//! Parses `.cf_seg` (unitig sequences), `.cf_seq` (reference tilings), and
//! `.json` (metadata including short references) to build the SSHash dictionary,
//! contig table, reference info, and optional equivalence class table.
//!
//! Corresponds to the C++ `build_contig_table()` in `build_contig_table.cpp`
//! plus dictionary construction from `build.cpp`.

use anyhow::{Context, Result, bail};
use sshash_lib::{BuildConfiguration, DictionaryBuilder, parse_cf_seg};
use std::collections::HashMap;
use std::io::BufRead;
use std::path::PathBuf;
use tracing::info;

use super::contig_table::ContigTableBuilder;
use super::eq_classes::{EqClassMapBuilder, Orientation};
use super::reference_index::ReferenceIndex;
use super::refinfo::RefInfo;

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------

/// Configuration for the index build pipeline.
pub struct BuildConfig {
    /// Basename for cuttlefish output files (.cf_seg, .cf_seq, .json).
    pub input_prefix: PathBuf,
    /// Output prefix for index files (.ssi, .ctab, .refinfo, .ectab).
    pub output_prefix: PathBuf,
    /// K-mer length.
    pub k: usize,
    /// Minimizer length.
    pub m: usize,
    /// Whether to build the equivalence class table.
    pub build_ec_table: bool,
    /// Number of threads (0 = all cores).
    pub num_threads: usize,
    /// Canonical k-mer mode.
    pub canonical: bool,
    /// Hash seed for the dictionary.
    pub seed: u64,
}

// ---------------------------------------------------------------------------
// Internal types
// ---------------------------------------------------------------------------

/// Segment metadata collected from .cf_seg parsing.
struct SegmentInfo {
    /// Index in file order (0-based rank).
    rank: u32,
    /// Sequence length in bases.
    len: u32,
}

/// Parsed tile from a .cf_seq tiling line.
enum Tile {
    /// A normal segment tile with ID and orientation.
    Segment { id: u64, is_fw: bool },
    /// A gap (N-tile) with the number of N's.
    Gap { count: u64 },
}

/// Parse a single tile token from a .cf_seq tiling line.
///
/// Formats: `<id>+` (forward), `<id>-` (reverse), `N<count>` (gap).
fn parse_tile(tok: &str) -> Result<Tile> {
    let last = tok.as_bytes().last()
        .context("empty tile token")?;

    match *last {
        b'+' => {
            let id: u64 = tok[..tok.len() - 1].parse()
                .with_context(|| format!("invalid segment id in '{tok}'"))?;
            Ok(Tile::Segment { id, is_fw: true })
        }
        b'-' => {
            let id: u64 = tok[..tok.len() - 1].parse()
                .with_context(|| format!("invalid segment id in '{tok}'"))?;
            Ok(Tile::Segment { id, is_fw: false })
        }
        _ => {
            if tok.starts_with('N') {
                let count: u64 = tok[1..].parse()
                    .with_context(|| format!("invalid N-tile count in '{tok}'"))?;
                Ok(Tile::Gap { count })
            } else {
                bail!("invalid tile token: '{tok}' (must end with +/- or start with N)")
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Main entry point
// ---------------------------------------------------------------------------

/// Build a piscem index from cuttlefish output files.
///
/// This is the main entry point for the build pipeline. It:
/// 1. Parses the `.cf_seg` file for unitig sequences
/// 2. Builds the SSHash dictionary
/// 3. Two-pass parses the `.cf_seq` file to build reference info and contig table
/// 4. Optionally builds the equivalence class table
/// 5. Assembles and saves the complete index
pub fn build_index(config: &BuildConfig) -> Result<()> {
    let k = config.k;
    let input = config.input_prefix.display().to_string();

    let cf_seg_path = format!("{input}.cf_seg");
    let cf_seq_path = format!("{input}.cf_seq");
    let json_path = format!("{input}.json");

    // Step 1: Parse .cf_seg — collect segment IDs, sequences, and lengths
    info!("Parsing segment file: {}", cf_seg_path);
    let cf_data = parse_cf_seg(&cf_seg_path)
        .with_context(|| format!("failed to parse {cf_seg_path}"))?;

    let num_segments = cf_data.len();
    info!("  {num_segments} segments");

    // Build segment info map: segment_id → SegmentInfo
    let mut id_to_info: HashMap<u64, SegmentInfo> = HashMap::with_capacity(num_segments);
    for (rank, (seg_id, seq)) in cf_data
        .segment_ids
        .iter()
        .zip(cf_data.sequences.iter())
        .enumerate()
    {
        id_to_info.insert(
            *seg_id,
            SegmentInfo {
                rank: rank as u32,
                len: seq.len() as u32,
            },
        );
    }

    // Step 2: Build SSHash dictionary from unitig sequences
    info!("Building SSHash dictionary (k={}, m={})", config.k, config.m);
    let mut build_cfg = BuildConfiguration::new(config.k, config.m)
        .map_err(|e| anyhow::anyhow!("invalid build configuration: {e}"))?;
    build_cfg.canonical = config.canonical;
    build_cfg.seed = config.seed;
    build_cfg.num_threads = config.num_threads;

    let dict_builder = DictionaryBuilder::new(build_cfg)
        .map_err(|e| anyhow::anyhow!("failed to create dictionary builder: {e}"))?;
    // Moves sequences into the builder
    let dict = dict_builder
        .build_from_sequences(cf_data.sequences)
        .map_err(|e| anyhow::anyhow!("failed to build dictionary: {e}"))?;
    info!("  Dictionary built: {} strings", dict.num_strings());

    // Step 3: First pass over .cf_seq — collect reference names, lengths
    info!("First pass over {cf_seq_path}");
    let (ref_names, ref_lens, max_ref_len) =
        first_pass_cf_seq(&cf_seq_path, &json_path, k, &id_to_info)?;
    let num_refs = ref_names.len();
    info!("  {num_refs} references, max_ref_len={max_ref_len}");

    // Step 4: Second pass over .cf_seq — build ContigTable
    info!("Second pass over {cf_seq_path} to fill contig entries");
    let contig_table = {
        let mut ctab_builder =
            ContigTableBuilder::new(num_segments, max_ref_len, num_refs as u64);
        second_pass_cf_seq(&cf_seq_path, k, &id_to_info, &mut ctab_builder)?;
        ctab_builder.build()
    };
    info!(
        "  ContigTable: {} contigs, {} entries",
        contig_table.num_contigs(),
        contig_table.num_entries()
    );

    // Step 5: Build EqClassMap (optional)
    let ec_table = if config.build_ec_table {
        info!("Building equivalence class table");
        let encoding = contig_table.encoding();
        let mut ec_builder = EqClassMapBuilder::new(num_segments);

        for tile_idx in 0..num_segments {
            let span = contig_table.contig_entries(tile_idx as u64);
            let mut label: Vec<(u32, Orientation)> = Vec::new();
            let mut prev_tid: u32 = 0;
            let mut prev_dir = Orientation::Forward;
            let mut first = true;

            for ce in span.iter() {
                let tid = encoding.transcript_id(ce);
                let dir = if encoding.orientation(ce) {
                    Orientation::Forward
                } else {
                    Orientation::ReverseComplement
                };

                if first || tid != prev_tid {
                    label.push((tid, dir));
                } else if dir != prev_dir {
                    // Same tid, different orientation → merge to BOTH
                    label.last_mut().unwrap().1 = Orientation::Both;
                }

                prev_tid = tid;
                prev_dir = dir;
                first = false;
            }

            ec_builder.add_tile(tile_idx, label);
        }

        let ec = ec_builder.build();
        info!(
            "  EqClassMap: {} tiles, {} ECs, {} label entries",
            ec.num_tiles(),
            ec.num_ecs(),
            ec.num_label_entries()
        );
        Some(ec)
    } else {
        None
    };

    // Step 6: Assemble and save
    let ref_info = RefInfo::new(ref_names, ref_lens);
    let index = ReferenceIndex::from_parts(dict, contig_table, ref_info, ec_table, None);

    // Create output directory if needed
    if let Some(parent) = config.output_prefix.parent() {
        if !parent.as_os_str().is_empty() && !parent.exists() {
            std::fs::create_dir_all(parent)
                .with_context(|| format!("failed to create output directory: {}", parent.display()))?;
        }
    }

    index.save(&config.output_prefix)?;
    info!("Index built and saved to {}", config.output_prefix.display());
    Ok(())
}

// ---------------------------------------------------------------------------
// First pass: reference info collection
// ---------------------------------------------------------------------------

/// First pass over the .cf_seq file: collects reference names and lengths.
///
/// Each line has format: `Reference:N_Sequence:<name>\t<tiles...>`
///
/// After processing .cf_seq, parses the .json file for short references
/// (sequences shorter than k that don't appear in the tiling).
fn first_pass_cf_seq(
    cf_seq_path: &str,
    json_path: &str,
    k: usize,
    id_to_info: &HashMap<u64, SegmentInfo>,
) -> Result<(Vec<String>, Vec<u64>, u64)> {
    let file = std::fs::File::open(cf_seq_path)
        .with_context(|| format!("failed to open {cf_seq_path}"))?;
    let reader = std::io::BufReader::new(file);

    let k_minus_1 = k as u64 - 1;
    let mut ref_names: Vec<String> = Vec::new();
    let mut ref_lens: Vec<u64> = Vec::new();
    let mut max_ref_len: u64 = 0;

    for (line_num, line_result) in reader.lines().enumerate() {
        let line = line_result
            .with_context(|| format!("failed to read line {} of {cf_seq_path}", line_num + 1))?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }

        // Split on tab: header \t tiles
        let (header, tiles_str) = line
            .split_once('\t')
            .with_context(|| format!("line {} missing tab separator", line_num + 1))?;

        // Extract reference name: everything after "Sequence:"
        let name = header
            .find("Sequence:")
            .map(|pos| &header[pos + 9..])
            .with_context(|| format!("line {} header missing 'Sequence:'", line_num + 1))?;

        // Parse tiles to compute reference length
        let mut current_offset: u64 = 0;
        for tok in tiles_str.split_whitespace() {
            match parse_tile(tok)? {
                Tile::Segment { id, .. } => {
                    let info = id_to_info.get(&id).with_context(|| {
                        format!("unknown segment id {id} on line {}", line_num + 1)
                    })?;
                    current_offset += info.len as u64 - k_minus_1;
                }
                Tile::Gap { count } => {
                    if current_offset > 0 {
                        current_offset += k_minus_1;
                    }
                    current_offset += count;
                }
            }
        }

        let ref_len = current_offset + k_minus_1;
        max_ref_len = max_ref_len.max(ref_len);
        ref_names.push(name.to_string());
        ref_lens.push(ref_len);

        if ref_names.len() % 10000 == 0 {
            info!("  processed {} references", ref_names.len());
        }
    }

    // Parse JSON for short references
    parse_short_refs(json_path, &mut ref_names, &mut ref_lens, &mut max_ref_len)?;

    Ok((ref_names, ref_lens, max_ref_len))
}

// ---------------------------------------------------------------------------
// Second pass: contig table population
// ---------------------------------------------------------------------------

/// Second pass over the .cf_seq file: populates the ContigTableBuilder with
/// occurrence records for each segment tile.
fn second_pass_cf_seq(
    cf_seq_path: &str,
    k: usize,
    id_to_info: &HashMap<u64, SegmentInfo>,
    builder: &mut ContigTableBuilder,
) -> Result<()> {
    let file = std::fs::File::open(cf_seq_path)
        .with_context(|| format!("failed to open {cf_seq_path}"))?;
    let reader = std::io::BufReader::new(file);

    let k_minus_1 = k as u64 - 1;
    let mut ref_id: u32 = 0;

    for line_result in reader.lines() {
        let line = line_result?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }

        // Split on tab — skip lines without tiles
        let tiles_str = match line.split_once('\t') {
            Some((_header, tiles)) => tiles,
            None => continue,
        };

        let mut current_offset: u64 = 0;
        for tok in tiles_str.split_whitespace() {
            match parse_tile(tok)? {
                Tile::Segment { id, is_fw } => {
                    let info = id_to_info
                        .get(&id)
                        .ok_or_else(|| anyhow::anyhow!("unknown segment id {id}"))?;
                    builder.add_occurrence(info.rank, ref_id, current_offset as u32, is_fw);
                    current_offset += info.len as u64 - k_minus_1;
                }
                Tile::Gap { count } => {
                    if current_offset > 0 {
                        current_offset += k_minus_1;
                    }
                    current_offset += count;
                }
            }
        }

        if ref_id % 10000 == 0 {
            info!("  processing reference #{ref_id}");
        }

        ref_id += 1;
    }

    Ok(())
}

// ---------------------------------------------------------------------------
// JSON short references
// ---------------------------------------------------------------------------

/// Parse short references from the cuttlefish .json metadata file.
///
/// Short references (shorter than k) are not part of the tiling and are
/// stored separately in the JSON under the `"short seqs"` key.
fn parse_short_refs(
    json_path: &str,
    ref_names: &mut Vec<String>,
    ref_lens: &mut Vec<u64>,
    max_ref_len: &mut u64,
) -> Result<()> {
    let file = match std::fs::File::open(json_path) {
        Ok(f) => f,
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => {
            info!("  No JSON metadata file found at {json_path}");
            return Ok(());
        }
        Err(e) => return Err(e).with_context(|| format!("failed to open {json_path}")),
    };

    let json: serde_json::Value = serde_json::from_reader(std::io::BufReader::new(file))
        .with_context(|| format!("failed to parse {json_path}"))?;

    // Check for "short seqs" (cuttlefish 2.x) or "short refs" (older versions)
    let short_refs = json.get("short seqs").or_else(|| json.get("short refs"));

    if let Some(entries) = short_refs.and_then(|v| v.as_array()) {
        for entry in entries {
            let arr = entry
                .as_array()
                .context("short seqs entry is not an array")?;
            if arr.len() != 2 {
                bail!("short seqs entry has {} elements, expected 2", arr.len());
            }
            let name = arr[0]
                .as_str()
                .context("short seqs name is not a string")?;
            let len = arr[1]
                .as_u64()
                .context("short seqs length is not a number")?;
            ref_names.push(name.to_string());
            ref_lens.push(len);
            *max_ref_len = (*max_ref_len).max(len);
        }
        info!("  {} short references from JSON", entries.len());
    }

    Ok(())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_tile_forward() {
        match parse_tile("42+").unwrap() {
            Tile::Segment { id, is_fw } => {
                assert_eq!(id, 42);
                assert!(is_fw);
            }
            _ => panic!("expected Segment"),
        }
    }

    #[test]
    fn test_parse_tile_reverse() {
        match parse_tile("137-").unwrap() {
            Tile::Segment { id, is_fw } => {
                assert_eq!(id, 137);
                assert!(!is_fw);
            }
            _ => panic!("expected Segment"),
        }
    }

    #[test]
    fn test_parse_tile_gap() {
        match parse_tile("N100").unwrap() {
            Tile::Gap { count } => assert_eq!(count, 100),
            _ => panic!("expected Gap"),
        }
    }

    #[test]
    fn test_parse_tile_large_id() {
        match parse_tile("93582850-").unwrap() {
            Tile::Segment { id, is_fw } => {
                assert_eq!(id, 93582850);
                assert!(!is_fw);
            }
            _ => panic!("expected Segment"),
        }
    }

    #[test]
    fn test_parse_tile_invalid() {
        assert!(parse_tile("xyz").is_err());
        assert!(parse_tile("").is_err());
    }

    #[test]
    fn test_offset_calculation() {
        // Simulate a reference with segments and an N-tile.
        // k = 5 for simplicity.
        let k_minus_1: u64 = 4;

        // Segment A (len=10): contribution = 10 - 4 = 6 → offset = 6
        // Segment B (len=8):  contribution =  8 - 4 = 4 → offset = 10
        // N5:                 offset > 0, so += 4 + 5 = 9 → offset = 19
        // Segment C (len=6):  contribution =  6 - 4 = 2 → offset = 21
        // ref_len = 21 + 4 = 25

        let mut current_offset: u64 = 0;

        // Segment A (len=10)
        current_offset += 10 - k_minus_1;
        assert_eq!(current_offset, 6);

        // Segment B (len=8)
        current_offset += 8 - k_minus_1;
        assert_eq!(current_offset, 10);

        // N5
        if current_offset > 0 {
            current_offset += k_minus_1;
        }
        current_offset += 5;
        assert_eq!(current_offset, 19);

        // Segment C (len=6)
        current_offset += 6 - k_minus_1;
        assert_eq!(current_offset, 21);

        let ref_len = current_offset + k_minus_1;
        assert_eq!(ref_len, 25);
    }

    #[test]
    fn test_offset_leading_n_tile() {
        // When a reference starts with an N-tile, we should NOT add k-1 overlap.
        let k_minus_1: u64 = 30;

        let mut current_offset: u64 = 0;

        // N10 at the start: current_offset == 0, so no overlap added
        if current_offset > 0 {
            current_offset += k_minus_1;
        }
        current_offset += 10;
        assert_eq!(current_offset, 10);

        // Segment (len=50): contribution = 50 - 30 = 20
        current_offset += 50 - k_minus_1;
        assert_eq!(current_offset, 30);

        let ref_len = current_offset + k_minus_1;
        assert_eq!(ref_len, 60);
    }

    #[test]
    fn test_short_refs_json() {
        let dir = std::env::temp_dir().join("piscem_rs_test_build");
        std::fs::create_dir_all(&dir).unwrap();
        let json_path = dir.join("test_short_refs.json");

        let json = r#"{
            "parameters info": {"k": 31},
            "short seqs": [
                ["short_ref_1", 12],
                ["short_ref_2", 25]
            ]
        }"#;
        std::fs::write(&json_path, json).unwrap();

        let mut names = vec!["existing".to_string()];
        let mut lens = vec![1000u64];
        let mut max_len = 1000u64;

        parse_short_refs(json_path.to_str().unwrap(), &mut names, &mut lens, &mut max_len)
            .unwrap();

        assert_eq!(names.len(), 3);
        assert_eq!(names[1], "short_ref_1");
        assert_eq!(names[2], "short_ref_2");
        assert_eq!(lens[1], 12);
        assert_eq!(lens[2], 25);
        // Short refs are smaller than existing max
        assert_eq!(max_len, 1000);

        // Cleanup
        std::fs::remove_file(&json_path).ok();
    }

    #[test]
    fn test_short_refs_json_missing_file() {
        let mut names = Vec::new();
        let mut lens = Vec::new();
        let mut max_len = 0u64;

        // Should not error if file doesn't exist
        parse_short_refs("/nonexistent/path.json", &mut names, &mut lens, &mut max_len).unwrap();

        assert!(names.is_empty());
        assert!(lens.is_empty());
    }

    #[test]
    fn test_short_refs_json_no_short_seqs() {
        let dir = std::env::temp_dir().join("piscem_rs_test_build");
        std::fs::create_dir_all(&dir).unwrap();
        let json_path = dir.join("test_no_short_refs.json");

        let json = r#"{"parameters info": {"k": 31}}"#;
        std::fs::write(&json_path, json).unwrap();

        let mut names = Vec::new();
        let mut lens = Vec::new();
        let mut max_len = 0u64;

        parse_short_refs(json_path.to_str().unwrap(), &mut names, &mut lens, &mut max_len)
            .unwrap();

        assert!(names.is_empty());

        std::fs::remove_file(&json_path).ok();
    }

    #[test]
    #[ignore] // Requires test_data directory and takes time to build dictionary
    fn test_build_from_test_data() {
        let output_dir = std::env::temp_dir().join("piscem_rs_build_integration_test");
        let config = BuildConfig {
            input_prefix: PathBuf::from(
                "test_data/gencode_pc_v44_dbg/gencode_pc_v44_index_cfish",
            ),
            output_prefix: output_dir.join("gencode_pc_v44_index"),
            k: 31,
            m: 19,
            build_ec_table: true,
            num_threads: 0,
            canonical: false,
            seed: 1,
        };

        build_index(&config).expect("build_index failed");

        // Verify we can load the index back
        let index = ReferenceIndex::load(&config.output_prefix, true, false)
            .expect("failed to reload index");

        assert!(index.num_contigs() > 0);
        assert!(index.num_refs() > 0);
        assert!(index.has_ec_table());

        // Cleanup
        std::fs::remove_dir_all(&output_dir).ok();
    }
}
