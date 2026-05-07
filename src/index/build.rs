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

use super::contig_table::ContigTableDirectBuilder;
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
    /// Use a single monolithic MPHF instead of partitioned.
    pub single_mphf: bool,
    /// Whether to emit TinyDictionary/TinyContigTable artifacts
    /// (`.tdct` + `.tct`). `Some(true)`/`Some(false)` forces the choice;
    /// `None` means "auto" — decide based on the canonical-k-mer count
    /// against [`crate::cli::AUTO_TINY_KMER_THRESHOLD`].
    pub emit_tiny: Option<bool>,
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
    /// Number of occurrences across all references (counted in first pass).
    count: u32,
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
    let last = tok.as_bytes().last().context("empty tile token")?;

    match *last {
        b'+' => {
            let id: u64 = tok[..tok.len() - 1]
                .parse()
                .with_context(|| format!("invalid segment id in '{tok}'"))?;
            Ok(Tile::Segment { id, is_fw: true })
        }
        b'-' => {
            let id: u64 = tok[..tok.len() - 1]
                .parse()
                .with_context(|| format!("invalid segment id in '{tok}'"))?;
            Ok(Tile::Segment { id, is_fw: false })
        }
        _ => {
            if let Some(rest) = tok.strip_prefix('N') {
                let count: u64 = rest
                    .parse()
                    .with_context(|| format!("invalid N-tile count in '{tok}'"))?;
                Ok(Tile::Gap { count })
            } else {
                bail!("invalid tile token: '{tok}' (must end with +/- or start with N)")
            }
        }
    }
}

#[derive(Default)]
struct ReferenceMetadata {
    short_refs: Vec<(String, u64)>,
    untiled_refs: HashMap<String, u64>,
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
    let mut cf_data =
        parse_cf_seg(&cf_seg_path).with_context(|| format!("failed to parse {cf_seg_path}"))?;

    let num_segments = cf_data.len();
    info!("  {num_segments} segments");

    // Extract segment IDs and lengths before moving sequences to dictionary builder.
    // We defer building the full id_to_info HashMap (~2-4 GB for 50M segments) until
    // AFTER dictionary build to reduce peak RSS during the memory-intensive MPHF phase.
    let segment_lens: Vec<u32> = cf_data.sequences.iter().map(|s| s.len() as u32).collect();
    let segment_ids = std::mem::take(&mut cf_data.segment_ids);

    // Step 2: Build SSHash dictionary from unitig sequences
    info!(
        "Building SSHash dictionary (k={}, m={})",
        config.k, config.m
    );
    let mut build_cfg = BuildConfiguration::new(config.k, config.m)
        .map_err(|e| anyhow::anyhow!("invalid build configuration: {e}"))?;
    build_cfg.canonical = config.canonical;
    build_cfg.seed = config.seed;
    build_cfg.num_threads = config.num_threads;
    build_cfg.partitioned_mphf = !config.single_mphf;

    let dict_builder = DictionaryBuilder::new(build_cfg)
        .map_err(|e| anyhow::anyhow!("failed to create dictionary builder: {e}"))?;
    // Moves sequences into the builder — cf_data.sequences consumed here
    let dict = dict_builder
        .build_from_sequences(cf_data.sequences)
        .map_err(|e| anyhow::anyhow!("failed to build dictionary: {e}"))?;
    info!("  Dictionary built: {} strings", dict.num_strings());

    // Now build id_to_info from the lightweight arrays (after peak RSS has passed)
    let mut id_to_info: HashMap<u64, SegmentInfo> = HashMap::with_capacity(num_segments);
    for (rank, (&seg_id, &len)) in segment_ids.iter().zip(segment_lens.iter()).enumerate() {
        id_to_info.insert(
            seg_id,
            SegmentInfo {
                rank: rank as u32,
                len,
                count: 0,
            },
        );
    }
    drop(segment_ids);
    drop(segment_lens);

    // Step 3: First pass over .cf_seq — collect reference names, lengths,
    //         and count per-unitig occurrences for the direct-write builder.
    info!("First pass over {cf_seq_path}");
    let (ref_names, ref_lens, max_ref_len) =
        first_pass_cf_seq(&cf_seq_path, &json_path, k, &mut id_to_info)?;
    let num_refs = ref_names.len();
    info!("  {num_refs} references, max_ref_len={max_ref_len}");

    // Step 4: Second pass over .cf_seq — build ContigTable using direct-write builder
    info!("Second pass over {cf_seq_path} to fill contig entries");
    let contig_table = {
        // Extract per-rank counts from id_to_info
        let mut counts = vec![0u32; num_segments];
        for info in id_to_info.values() {
            counts[info.rank as usize] = info.count;
        }
        let mut ctab_builder = ContigTableDirectBuilder::new(&counts, max_ref_len, num_refs as u64);
        second_pass_cf_seq(&cf_seq_path, k, &id_to_info, &mut ctab_builder)?;
        ctab_builder.build()
    };
    info!(
        "  ContigTable: {} contigs, {} entries",
        contig_table.num_contigs(),
        contig_table.num_entries()
    );
    // Free segment metadata — no longer needed after both passes
    drop(id_to_info);

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
    if let Some(parent) = config.output_prefix.parent()
        && !parent.as_os_str().is_empty()
        && !parent.exists()
    {
        std::fs::create_dir_all(parent)
            .with_context(|| format!("failed to create output directory: {}", parent.display()))?;
    }

    index.save(&config.output_prefix)?;
    info!(
        "Index built and saved to {}",
        config.output_prefix.display()
    );

    let emit_tiny = match config.emit_tiny {
        Some(v) => v,
        None => {
            let dict = index.dict();
            let total_bases = dict.spss().total_bases();
            let num_kmers = total_bases.saturating_sub(dict.num_strings() * (k as u64 - 1));
            let below = num_kmers < crate::cli::AUTO_TINY_KMER_THRESHOLD;
            info!(
                "auto dict selection: {} canonical k-mers ({} threshold), \
                 emitting Tiny artifacts: {}",
                num_kmers,
                crate::cli::AUTO_TINY_KMER_THRESHOLD,
                below
            );
            below
        }
    };
    if emit_tiny {
        use sshash_lib::dispatch_on_k;
        info!("Emitting Tiny index artifacts (.tdct + .tct)");
        dispatch_on_k!(k, K => {
            index.save_tiny_artifacts::<K>(&config.output_prefix)?;
        });
    }
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
    id_to_info: &mut HashMap<u64, SegmentInfo>,
) -> Result<(Vec<String>, Vec<u64>, u64)> {
    let metadata = parse_reference_metadata(json_path)?;
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
        let line = line.trim_end_matches(['\n', '\r']);
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

        // Parse tiles to compute reference length and count per-unitig occurrences
        let mut current_offset: u64 = 0;
        let mut saw_tile = false;
        for tok in tiles_str.split_whitespace() {
            saw_tile = true;
            match parse_tile(tok)? {
                Tile::Segment { id, .. } => {
                    let info = id_to_info.get_mut(&id).with_context(|| {
                        format!("unknown segment id {id} on line {}", line_num + 1)
                    })?;
                    info.count += 1;
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

        let ref_len = if saw_tile {
            current_offset + k_minus_1
        } else {
            metadata.untiled_refs.get(name).copied().with_context(|| {
                format!(
                    "line {} has no tiles; expected '{}'-length entry in JSON metadata key 'untiled seqs'",
                    line_num + 1,
                    name
                )
            })?
        };
        max_ref_len = max_ref_len.max(ref_len);
        ref_names.push(name.to_string());
        ref_lens.push(ref_len);

        if ref_names.len().is_multiple_of(10000) {
            info!("  processed {} references", ref_names.len());
        }
    }

    // Append short references from JSON metadata.
    append_short_refs(&metadata, &mut ref_names, &mut ref_lens, &mut max_ref_len);

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
    builder: &mut ContigTableDirectBuilder,
) -> Result<()> {
    let file = std::fs::File::open(cf_seq_path)
        .with_context(|| format!("failed to open {cf_seq_path}"))?;
    let reader = std::io::BufReader::new(file);

    let k_minus_1 = k as u64 - 1;
    let mut ref_id: u32 = 0;

    for line_result in reader.lines() {
        let line = line_result?;
        let line = line.trim_end_matches(['\n', '\r']);
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

        if ref_id.is_multiple_of(10000) {
            info!("  processing reference #{ref_id}");
        }

        ref_id += 1;
    }

    Ok(())
}

// ---------------------------------------------------------------------------
// JSON short references
// ---------------------------------------------------------------------------

/// Parse reference metadata from the cuttlefish .json sidecar.
///
/// Supported keys:
/// - `"short seqs"` / `"short refs"`: references shorter than `k`
/// - `"untiled seqs"`: references present in `.cf_seq` but with no tile list
fn parse_reference_metadata(json_path: &str) -> Result<ReferenceMetadata> {
    let file = match std::fs::File::open(json_path) {
        Ok(f) => f,
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => {
            info!("  No JSON metadata file found at {json_path}");
            return Ok(ReferenceMetadata::default());
        }
        Err(e) => return Err(e).with_context(|| format!("failed to open {json_path}")),
    };

    let json: serde_json::Value = serde_json::from_reader(std::io::BufReader::new(file))
        .with_context(|| format!("failed to parse {json_path}"))?;

    let mut metadata = ReferenceMetadata::default();

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
            let name = arr[0].as_str().context("short seqs name is not a string")?;
            let len = arr[1]
                .as_u64()
                .context("short seqs length is not a number")?;
            metadata.short_refs.push((name.to_string(), len));
        }
        info!("  {} short references from JSON", entries.len());
    }

    if let Some(entries) = json.get("untiled seqs").and_then(|v| v.as_array()) {
        for entry in entries {
            let arr = entry
                .as_array()
                .context("untiled seqs entry is not an array")?;
            if arr.len() != 2 {
                bail!("untiled seqs entry has {} elements, expected 2", arr.len());
            }
            let name = arr[0]
                .as_str()
                .context("untiled seqs name is not a string")?;
            let len = arr[1]
                .as_u64()
                .context("untiled seqs length is not a number")?;
            metadata.untiled_refs.insert(name.to_string(), len);
        }
        info!("  {} untiled references from JSON", entries.len());
    }

    Ok(metadata)
}

fn append_short_refs(
    metadata: &ReferenceMetadata,
    ref_names: &mut Vec<String>,
    ref_lens: &mut Vec<u64>,
    max_ref_len: &mut u64,
) {
    for (name, len) in &metadata.short_refs {
        ref_names.push(name.clone());
        ref_lens.push(*len);
        *max_ref_len = (*max_ref_len).max(*len);
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

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
    fn test_reference_metadata_json() {
        let dir = std::env::temp_dir().join("piscem_rs_test_build");
        std::fs::create_dir_all(&dir).unwrap();
        let json_path = dir.join("test_reference_metadata.json");

        let json = r#"{
            "parameters info": {"k": 31},
            "short seqs": [
                ["short_ref_1", 12],
                ["short_ref_2", 25]
            ],
            "untiled seqs": [
                ["all_n_1", 337],
                ["all_n_2", 91]
            ]
        }"#;
        std::fs::write(&json_path, json).unwrap();

        let metadata = parse_reference_metadata(json_path.to_str().unwrap()).unwrap();

        assert_eq!(metadata.short_refs.len(), 2);
        assert_eq!(metadata.short_refs[0], ("short_ref_1".to_string(), 12));
        assert_eq!(metadata.short_refs[1], ("short_ref_2".to_string(), 25));
        assert_eq!(metadata.untiled_refs.get("all_n_1"), Some(&337));
        assert_eq!(metadata.untiled_refs.get("all_n_2"), Some(&91));

        // Cleanup
        std::fs::remove_file(&json_path).ok();
    }

    #[test]
    fn test_reference_metadata_json_missing_file() {
        let metadata = parse_reference_metadata("/nonexistent/path.json").unwrap();
        assert!(metadata.short_refs.is_empty());
        assert!(metadata.untiled_refs.is_empty());
    }

    #[test]
    fn test_reference_metadata_json_no_short_seqs() {
        let dir = std::env::temp_dir().join("piscem_rs_test_build");
        std::fs::create_dir_all(&dir).unwrap();
        let json_path = dir.join("test_no_short_refs.json");

        let json = r#"{"parameters info": {"k": 31}}"#;
        std::fs::write(&json_path, json).unwrap();

        let metadata = parse_reference_metadata(json_path.to_str().unwrap()).unwrap();
        assert!(metadata.short_refs.is_empty());
        assert!(metadata.untiled_refs.is_empty());

        std::fs::remove_file(&json_path).ok();
    }

    #[test]
    fn test_first_pass_cf_seq_preserves_empty_tiles_field() {
        let dir = tempdir().unwrap();
        let cf_seq_path = dir.path().join("empty_tiles.cf_seq");
        let json_path = dir.path().join("empty_tiles.json");

        std::fs::write(&cf_seq_path, "Reference:N_Sequence:all_n\t\n").unwrap();
        std::fs::write(&json_path, r#"{"untiled seqs":[["all_n",337]]}"#).unwrap();

        let mut id_to_info = HashMap::new();
        let (names, lens, max_ref_len) = first_pass_cf_seq(
            cf_seq_path.to_str().unwrap(),
            json_path.to_str().unwrap(),
            31,
            &mut id_to_info,
        )
        .unwrap();

        assert_eq!(names, vec!["all_n".to_string()]);
        assert_eq!(lens, vec![337]);
        assert_eq!(max_ref_len, 337);
    }

    #[test]
    fn test_first_pass_cf_seq_rejects_empty_tiles_without_metadata() {
        let dir = tempdir().unwrap();
        let cf_seq_path = dir.path().join("missing_untiled_metadata.cf_seq");
        let json_path = dir.path().join("missing_untiled_metadata.json");

        std::fs::write(&cf_seq_path, "Reference:N_Sequence:all_n\t\n").unwrap();
        std::fs::write(&json_path, "{}").unwrap();

        let mut id_to_info = HashMap::new();
        let err = first_pass_cf_seq(
            cf_seq_path.to_str().unwrap(),
            json_path.to_str().unwrap(),
            31,
            &mut id_to_info,
        )
        .unwrap_err();

        let msg = format!("{err:#}");
        assert!(msg.contains("untiled seqs"));
    }

    #[test]
    fn test_second_pass_cf_seq_counts_empty_tiles_reference() {
        let dir = tempdir().unwrap();
        let cf_seq_path = dir.path().join("empty_tiles_then_segment.cf_seq");

        std::fs::write(
            &cf_seq_path,
            "Reference:N_Sequence:all_n\t\nReference:N_Sequence:segmented\t7+\n",
        )
        .unwrap();

        let mut id_to_info = HashMap::new();
        id_to_info.insert(
            7,
            SegmentInfo {
                rank: 0,
                len: 35,
                count: 1,
            },
        );

        let mut builder = ContigTableDirectBuilder::new(&[1], 64, 2);
        second_pass_cf_seq(cf_seq_path.to_str().unwrap(), 31, &id_to_info, &mut builder).unwrap();
        let contig_table = builder.build();
        let encoding = contig_table.encoding();
        let entries: Vec<u64> = contig_table.contig_entries(0).iter().collect();

        assert_eq!(entries.len(), 1);
        assert_eq!(encoding.transcript_id(entries[0]), 1);
        assert_eq!(encoding.pos(entries[0]), 0);
        assert!(encoding.orientation(entries[0]));
    }

    #[test]
    #[ignore] // Requires test_data directory and takes time to build dictionary
    fn test_build_from_test_data() {
        let output_dir = std::env::temp_dir().join("piscem_rs_build_integration_test");
        let config = BuildConfig {
            input_prefix: PathBuf::from("test_data/gencode_pc_v44_dbg/gencode_pc_v44_index_cfish"),
            output_prefix: output_dir.join("gencode_pc_v44_index"),
            k: 31,
            m: 19,
            build_ec_table: true,
            num_threads: 0,
            canonical: false,
            seed: 1,
            single_mphf: false,
            emit_tiny: Some(false),
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
