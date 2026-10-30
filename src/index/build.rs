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

use crate::hash::FixedMap;
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
    /// Output prefix for index files: each file is the prefix plus an
    /// extension (`<prefix>.ssi`, `<prefix>.ctab`, ...; see
    /// [`crate::index::paths`]). It must end in a file name
    /// ([`crate::index::paths::validate_index_prefix`]).
    pub output_prefix: PathBuf,
    /// K-mer length.
    pub k: usize,
    /// Minimizer length.
    pub m: usize,
    /// Whether to build the equivalence class table.
    pub build_ec_table: bool,
    /// Number of threads (0 = all cores).
    pub num_threads: usize,
    /// Hash seed for the dictionary.
    pub seed: u64,
    /// How sshash builds the minimizer MPHF (see [`sshash_lib::MphfBuild`]).
    pub mphf_build: sshash_lib::MphfBuild,
    /// Whether to emit TinyDictionary/TinyContigTable artifacts
    /// (`.tdct` + `.tct`). `Some(true)`/`Some(false)` forces the choice;
    /// `None` means "auto" — decide based on the canonical-k-mer count
    /// against [`crate::cli::AUTO_TINY_KMER_THRESHOLD`].
    pub emit_tiny: Option<bool>,
    /// Directory sshash uses for its temporary files (external minimizer
    /// sort, offline MPHF construction records).
    /// `None` keeps sshash's default (`sshash_tmp`, relative to the current
    /// working directory); `Some(dir)` redirects it so the caller controls the
    /// location and cleanup. The directory is created by sshash if absent.
    pub tmp_dir: Option<PathBuf>,
    /// RAM ceiling, in GiB, for sshash's external minimizer sort (the dominant
    /// build-time memory/disk trade-off).
    ///
    /// `Some(n)` is used verbatim, `Some(0)` meaning sshash's "unlimited".
    /// `None` auto-detects from the tightest memory limit the platform exposes
    /// — see [`crate::index::mem_budget`], which is careful to respect cgroup
    /// limits rather than trusting physical RAM on a scheduler-managed node.
    pub ram_limit_gib: Option<usize>,
}

// ---------------------------------------------------------------------------
// Internal types
// ---------------------------------------------------------------------------

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
    // Before any work: a prefix naming a directory would otherwise surface
    // only when the files are written, beside that directory.
    crate::index::paths::validate_index_prefix(&config.output_prefix)?;
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

    // The build's own data-parallel steps run on `num_threads` workers, like
    // the dictionary build (0 = all cores).
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(config.num_threads)
        .build()
        .context("failed to create the index-build thread pool")?;

    // Step 2: Parse .cf_seq — reference names and lengths, each reference's
    //         tiles, and per-unitig occurrence counts for the direct-write
    //         builder — with segments provisionally ranked in file order.
    info!("Parsing tilings in {cf_seq_path}");
    let mut tilings = {
        let mut id_to_rank: FixedMap<u64, u32> =
            FixedMap::with_capacity_and_hasher(num_segments, crate::hash::fixed_state());
        for (rank, &seg_id) in cf_data.segment_ids.iter().enumerate() {
            id_to_rank.insert(seg_id, rank as u32);
        }
        let file_lens: Vec<u32> = cf_data.sequences.iter().map(|s| s.len() as u32).collect();
        pool.install(|| parse_cf_seq(&cf_seq_path, &json_path, k, &id_to_rank, &file_lens))?
    };
    let num_refs = tilings.names.len();
    let max_ref_len = tilings.max_ref_len;
    info!("  {num_refs} references, max_ref_len={max_ref_len}");

    // Number the contigs in order of first appearance in the tilings, rather
    // than in .cf_seg file order: cDBG builders write segments in whatever
    // order their threads finish them (cf1-rs does), so file order made every
    // index file differ from build to build of the same reference. Order of
    // first appearance is reproducible and keeps the contigs of a reference
    // close together.
    let order = first_appearance_order(&tilings, &cf_data.segment_ids);
    pool.install(|| tilings.renumber(&order));
    let sequences: Vec<String> = {
        let mut by_file_rank = std::mem::take(&mut cf_data.sequences);
        order
            .iter()
            .map(|&old| std::mem::take(&mut by_file_rank[old as usize]))
            .collect()
    };
    drop(cf_data);
    let segment_lens: Vec<u32> = sequences.iter().map(|s| s.len() as u32).collect();

    // Step 3: Build SSHash dictionary from unitig sequences
    info!(
        "Building SSHash dictionary (k={}, m={})",
        config.k, config.m
    );
    let mut build_cfg = BuildConfiguration::new(config.k, config.m)
        .map_err(|e| anyhow::anyhow!("invalid build configuration: {e}"))?;
    build_cfg.seed = config.seed;
    build_cfg.num_threads = config.num_threads;
    build_cfg.mphf_build = config.mphf_build;
    if let Some(dir) = &config.tmp_dir {
        build_cfg.tmp_dirname = dir.clone();
    }
    build_cfg.ram_limit_gib = crate::index::mem_budget::resolve_ram_limit_gib(config.ram_limit_gib);

    let dict_builder = DictionaryBuilder::new(build_cfg)
        .map_err(|e| anyhow::anyhow!("failed to create dictionary builder: {e}"))?;
    // Moves sequences into the builder
    let dict = dict_builder
        .build_from_sequences(sequences)
        .map_err(|e| anyhow::anyhow!("failed to build dictionary: {e}"))?;
    info!("  Dictionary built: {} strings", dict.num_strings());

    // Step 4: Fill the ContigTable from the parsed tiles.
    info!("Filling contig entries");
    let contig_table = {
        let mut ctab_builder =
            ContigTableDirectBuilder::new(&tilings.counts, max_ref_len, num_refs as u64);
        fill_contig_table(&tilings, k, &segment_lens, &mut ctab_builder);
        ctab_builder.build()
    };
    info!(
        "  ContigTable: {} contigs, {} entries",
        contig_table.num_contigs(),
        contig_table.num_entries()
    );
    drop(segment_lens);
    let Tilings {
        names: ref_names,
        ref_lens,
        ..
    } = tilings;

    // Step 5: Build EqClassMap (optional)
    let ec_table = if config.build_ec_table {
        info!("Building equivalence class table");
        use rayon::prelude::*;

        let encoding = contig_table.encoding();
        let mut ec_builder = EqClassMapBuilder::new(num_segments);

        // A tile's sorted label (and its hash) depends only on its contig
        // entries, so blocks of labels are computed in parallel; they are then
        // interned in tile order, which keeps EC ids in first-occurrence order.
        let tile_label = |tile_idx: usize| {
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
            label.sort();
            label
        };
        const LABEL_BLOCK: usize = 1 << 18;
        for block_start in (0..num_segments).step_by(LABEL_BLOCK) {
            let block = block_start..(block_start + LABEL_BLOCK).min(num_segments);
            let labelled: Vec<_> = pool.install(|| {
                block
                    .clone()
                    .into_par_iter()
                    .map(|tile_idx| {
                        let label = tile_label(tile_idx);
                        let hash = ec_builder.label_hash(&label);
                        (label, hash)
                    })
                    .collect()
            });
            for (tile_idx, (label, hash)) in block.zip(labelled) {
                ec_builder.add_sorted_tile(tile_idx, label, hash);
            }
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
// Reference tilings (.cf_seq)
// ---------------------------------------------------------------------------

/// Marks a gap (N-run) in a packed tile code; the low bits hold its length.
const GAP_TILE: u64 = 1 << 63;

/// The reference tilings of a `.cf_seq` file, parsed once.
///
/// Tiles are packed one per `u64` and the per-reference lists concatenated
/// per parse chunk, so filling the contig table afterwards re-reads neither the
/// text nor the segment-id map.
struct Tilings {
    /// Reference names in file order, then the short references from the JSON.
    names: Vec<String>,
    /// Reference lengths, parallel to `names`.
    ref_lens: Vec<u64>,
    max_ref_len: u64,
    /// Occurrences of each unitig (by rank) across all tilings.
    counts: Vec<u32>,
    /// The tilings, in reference order, split into the chunks they were
    /// parsed in.
    chunks: Vec<TileChunk>,
}

/// The packed tiles of a run of consecutive references.
#[derive(Default)]
struct TileChunk {
    /// Segment tiles are `rank << 1 | is_fw`; gaps are `GAP_TILE | count`.
    tiles: Vec<u64>,
    /// `tiles[starts[i]..starts[i + 1]]` belong to the chunk's `i`-th reference.
    starts: Vec<usize>,
    names: Vec<String>,
    ref_lens: Vec<u64>,
}

/// Parse a decimal `u64` from ASCII digits, rejecting anything else.
fn parse_u64(digits: &[u8]) -> Option<u64> {
    if digits.is_empty() {
        return None;
    }
    digits.iter().try_fold(0u64, |acc, &b| {
        let d = b.wrapping_sub(b'0');
        (d < 10)
            .then_some(())
            .and_then(|()| acc.checked_mul(10)?.checked_add(d as u64))
    })
}

/// Parse the `.cf_seq` file: reference names and lengths, the packed tiles of
/// every reference, and how often each unitig occurs.
///
/// Each line has format: `Reference:N_Sequence:<name>\t<tiles...>`. Lines are
/// parsed in parallel in chunks; references keep their file order. Short
/// references (sequences shorter than `k`, absent from the tiling) are appended
/// from the `.json` metadata, and references with an empty tile list take their
/// length from it.
fn parse_cf_seq(
    cf_seq_path: &str,
    json_path: &str,
    k: usize,
    id_to_rank: &FixedMap<u64, u32>,
    segment_lens: &[u32],
) -> Result<Tilings> {
    use rayon::prelude::*;
    use std::sync::atomic::{AtomicU32, Ordering};

    let metadata = parse_reference_metadata(json_path)?;
    let text =
        std::fs::read(cf_seq_path).with_context(|| format!("failed to read {cf_seq_path}"))?;

    // Chunk boundaries at line starts, a few chunks per thread for balance.
    let target = text.len() / (rayon::current_num_threads() * 8).max(1) + 1;
    let mut bounds = vec![0usize];
    while let Some(&last) = bounds.last() {
        if last >= text.len() {
            break;
        }
        let probe = (last + target).min(text.len());
        let next = text[probe..]
            .iter()
            .position(|&b| b == b'\n')
            .map_or(text.len(), |p| probe + p + 1);
        bounds.push(next);
    }

    let counts: Vec<AtomicU32> = (0..segment_lens.len()).map(|_| AtomicU32::new(0)).collect();
    let k_minus_1 = k as u64 - 1;
    let chunks: Vec<TileChunk> = bounds
        .par_windows(2)
        .map(|w| {
            let mut chunk = TileChunk::default();
            for line in text[w[0]..w[1]].split(|&b| b == b'\n') {
                let line = line.strip_suffix(b"\r").unwrap_or(line);
                if line.is_empty() {
                    continue;
                }
                let line = std::str::from_utf8(line).context(".cf_seq line is not UTF-8")?;
                // Split on tab: header \t tiles
                let (header, tiles_str) = line
                    .split_once('\t')
                    .with_context(|| format!("line missing tab separator: '{line}'"))?;
                // Extract reference name: everything after "Sequence:"
                let name = header
                    .find("Sequence:")
                    .map(|pos| &header[pos + 9..])
                    .with_context(|| format!("header missing 'Sequence:': '{header}'"))?;

                chunk.starts.push(chunk.tiles.len());
                let mut current_offset: u64 = 0;
                let mut saw_tile = false;
                for tok in tiles_str.as_bytes().split(|&b| b == b' ') {
                    let Some(&last) = tok.last() else { continue };
                    saw_tile = true;
                    let code = match last {
                        b'+' | b'-' => {
                            let id = parse_u64(&tok[..tok.len() - 1]).with_context(|| {
                                format!("invalid segment id in '{}'", String::from_utf8_lossy(tok))
                            })?;
                            let rank = *id_to_rank.get(&id).with_context(|| {
                                format!("unknown segment id {id} in the tiling of '{name}'")
                            })?;
                            counts[rank as usize].fetch_add(1, Ordering::Relaxed);
                            current_offset += segment_lens[rank as usize] as u64 - k_minus_1;
                            (rank as u64) << 1 | (last == b'+') as u64
                        }
                        b'0'..=b'9' if tok[0] == b'N' => {
                            let count = parse_u64(&tok[1..]).with_context(|| {
                                format!(
                                    "invalid N-tile count in '{}'",
                                    String::from_utf8_lossy(tok)
                                )
                            })?;
                            if current_offset > 0 {
                                current_offset += k_minus_1;
                            }
                            current_offset += count;
                            GAP_TILE | count
                        }
                        _ => bail!(
                            "invalid tile token: '{}' (must end with +/- or start with N)",
                            String::from_utf8_lossy(tok)
                        ),
                    };
                    chunk.tiles.push(code);
                }

                let ref_len = if saw_tile {
                    current_offset + k_minus_1
                } else {
                    metadata.untiled_refs.get(name).copied().with_context(|| {
                        format!(
                            "reference '{name}' has no tiles; expected its length in JSON \
                             metadata key 'untiled seqs'"
                        )
                    })?
                };
                chunk.names.push(name.to_string());
                chunk.ref_lens.push(ref_len);
            }
            chunk.starts.push(chunk.tiles.len());
            Ok(chunk)
        })
        .collect::<Result<_>>()?;
    drop(text);

    let num_refs: usize = chunks.iter().map(|c| c.names.len()).sum();
    let mut names = Vec::with_capacity(num_refs + metadata.short_refs.len());
    let mut ref_lens = Vec::with_capacity(names.capacity());
    for chunk in &chunks {
        names.extend_from_slice(&chunk.names);
        ref_lens.extend_from_slice(&chunk.ref_lens);
    }
    let mut max_ref_len = ref_lens.iter().copied().max().unwrap_or(0);

    // Append short references from JSON metadata.
    append_short_refs(&metadata, &mut names, &mut ref_lens, &mut max_ref_len);

    Ok(Tilings {
        names,
        ref_lens,
        max_ref_len,
        counts: counts.into_iter().map(AtomicU32::into_inner).collect(),
        chunks,
    })
}

/// The order in which segments first appear in `tilings` (reference order,
/// then position), as a list of their current ranks. Segments that appear in
/// no tiling follow, by segment id, so the order never depends on the order
/// segments were listed in.
fn first_appearance_order(tilings: &Tilings, segment_ids: &[u64]) -> Vec<u32> {
    let mut seen = vec![false; segment_ids.len()];
    let mut order = Vec::with_capacity(segment_ids.len());
    for code in tilings.chunks.iter().flat_map(|c| &c.tiles) {
        if code & GAP_TILE == 0 {
            let rank = (code >> 1) as u32;
            if !std::mem::replace(&mut seen[rank as usize], true) {
                order.push(rank);
            }
        }
    }
    if order.len() < segment_ids.len() {
        let mut unseen: Vec<u32> = (0..segment_ids.len() as u32)
            .filter(|&r| !seen[r as usize])
            .collect();
        unseen.sort_unstable_by_key(|&r| segment_ids[r as usize]);
        order.extend(unseen);
    }
    order
}

impl Tilings {
    /// Re-rank the segments: `order[new_rank]` is the segment's current rank.
    fn renumber(&mut self, order: &[u32]) {
        use rayon::prelude::*;

        let mut new_rank = vec![0u32; order.len()];
        for (new, &old) in order.iter().enumerate() {
            new_rank[old as usize] = new as u32;
        }
        self.counts = order.iter().map(|&old| self.counts[old as usize]).collect();
        self.chunks.par_iter_mut().for_each(|chunk| {
            for code in &mut chunk.tiles {
                if *code & GAP_TILE == 0 {
                    *code = (new_rank[(*code >> 1) as usize] as u64) << 1 | (*code & 1);
                }
            }
        });
    }
}

/// Populate the ContigTableBuilder with one occurrence record per segment tile,
/// in reference order.
fn fill_contig_table(
    tilings: &Tilings,
    k: usize,
    segment_lens: &[u32],
    builder: &mut ContigTableDirectBuilder,
) {
    let k_minus_1 = k as u64 - 1;
    let mut ref_id: u32 = 0;
    for chunk in &tilings.chunks {
        for w in chunk.starts.windows(2) {
            let mut current_offset: u64 = 0;
            for &code in &chunk.tiles[w[0]..w[1]] {
                if code & GAP_TILE != 0 {
                    if current_offset > 0 {
                        current_offset += k_minus_1;
                    }
                    current_offset += code & !GAP_TILE;
                } else {
                    let rank = (code >> 1) as u32;
                    builder.add_occurrence(rank, ref_id, current_offset as u32, code & 1 == 1);
                    current_offset += segment_lens[rank as usize] as u64 - k_minus_1;
                }
            }
            ref_id += 1;
        }
    }
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
    fn test_parse_u64() {
        assert_eq!(parse_u64(b"0"), Some(0));
        assert_eq!(parse_u64(b"93582850"), Some(93582850));
        assert_eq!(parse_u64(b"18446744073709551615"), Some(u64::MAX));
        assert_eq!(parse_u64(b"18446744073709551616"), None);
        assert_eq!(parse_u64(b""), None);
        assert_eq!(parse_u64(b"12a"), None);
    }

    #[test]
    fn test_parse_cf_seq_preserves_empty_tiles_field() {
        let dir = tempdir().unwrap();
        let cf_seq_path = dir.path().join("empty_tiles.cf_seq");
        let json_path = dir.path().join("empty_tiles.json");

        std::fs::write(&cf_seq_path, "Reference:N_Sequence:all_n\t\n").unwrap();
        std::fs::write(&json_path, r#"{"untiled seqs":[["all_n",337]]}"#).unwrap();

        let id_to_rank = FixedMap::with_hasher(crate::hash::fixed_state());
        let tilings = parse_cf_seq(
            cf_seq_path.to_str().unwrap(),
            json_path.to_str().unwrap(),
            31,
            &id_to_rank,
            &[],
        )
        .unwrap();

        assert_eq!(tilings.names, vec!["all_n".to_string()]);
        assert_eq!(tilings.ref_lens, vec![337]);
        assert_eq!(tilings.max_ref_len, 337);
    }

    #[test]
    fn test_parse_cf_seq_rejects_empty_tiles_without_metadata() {
        let dir = tempdir().unwrap();
        let cf_seq_path = dir.path().join("missing_untiled_metadata.cf_seq");
        let json_path = dir.path().join("missing_untiled_metadata.json");

        std::fs::write(&cf_seq_path, "Reference:N_Sequence:all_n\t\n").unwrap();
        std::fs::write(&json_path, "{}").unwrap();

        let id_to_rank = FixedMap::with_hasher(crate::hash::fixed_state());
        let err = parse_cf_seq(
            cf_seq_path.to_str().unwrap(),
            json_path.to_str().unwrap(),
            31,
            &id_to_rank,
            &[],
        )
        .err()
        .unwrap();

        let msg = format!("{err:#}");
        assert!(msg.contains("untiled seqs"));
    }

    #[test]
    fn test_parse_cf_seq_rejects_bad_tiles() {
        let dir = tempdir().unwrap();
        let json_path = dir.path().join("none.json");
        let mut id_to_rank = FixedMap::with_hasher(crate::hash::fixed_state());
        id_to_rank.insert(7, 0);
        for (case, line) in [
            ("unknown", "Reference:1_Sequence:r\t8+\n"),
            ("garbage", "Reference:1_Sequence:r\txyz\n"),
            ("bad_gap", "Reference:1_Sequence:r\t7+ N1x\n"),
            ("no_tab", "Reference:1_Sequence:r 7+\n"),
        ] {
            let cf_seq_path = dir.path().join(format!("{case}.cf_seq"));
            std::fs::write(&cf_seq_path, line).unwrap();
            let res = parse_cf_seq(
                cf_seq_path.to_str().unwrap(),
                json_path.to_str().unwrap(),
                31,
                &id_to_rank,
                &[35],
            );
            assert!(res.is_err(), "{case} should be rejected");
        }
    }

    #[test]
    fn test_parse_and_fill_cf_seq_offsets_and_order() {
        // k = 5: segment 7 (rank 0, len 10), segment 9 (rank 1, len 8).
        let dir = tempdir().unwrap();
        let cf_seq_path = dir.path().join("tiles.cf_seq");
        let json_path = dir.path().join("tiles.json");
        std::fs::write(
            &cf_seq_path,
            "Reference:1_Sequence:all_n\t\n\
             Reference:1_Sequence:a\t7+ 9- N5 7-\r\n\
             \n\
             Reference:1_Sequence:b\tN3 9+\n",
        )
        .unwrap();
        std::fs::write(
            &json_path,
            r#"{"untiled seqs":[["all_n",40]],"short seqs":[["s",3]]}"#,
        )
        .unwrap();
        let mut id_to_rank = FixedMap::with_hasher(crate::hash::fixed_state());
        id_to_rank.insert(7, 0);
        id_to_rank.insert(9, 1);
        let lens = [10u32, 8];

        let tilings = parse_cf_seq(
            cf_seq_path.to_str().unwrap(),
            json_path.to_str().unwrap(),
            5,
            &id_to_rank,
            &lens,
        )
        .unwrap();
        assert_eq!(tilings.names, ["all_n", "a", "b", "s"]);
        // a: 6 + 4, N5 -> +4 +5, then 7- -> +6 => 25, + k-1 = 29
        // b: N3 -> 3, 9+ -> +4 => 7, + k-1 = 11
        assert_eq!(tilings.ref_lens, [40, 29, 11, 3]);
        assert_eq!(tilings.max_ref_len, 40);
        assert_eq!(tilings.counts, [2, 2]);

        let mut builder = ContigTableDirectBuilder::new(&tilings.counts, 40, 4);
        fill_contig_table(&tilings, 5, &lens, &mut builder);
        let table = builder.build();
        let enc = table.encoding();
        let decode = |rank: u64| -> Vec<(u32, u64, bool)> {
            table
                .contig_entries(rank)
                .iter()
                .map(|e| (enc.transcript_id(e), enc.pos(e) as u64, enc.orientation(e)))
                .collect()
        };
        assert_eq!(decode(0), [(1, 0, true), (1, 19, false)]);
        assert_eq!(decode(1), [(1, 6, false), (2, 3, true)]);
    }

    #[test]
    fn test_first_appearance_order_and_renumber() {
        let dir = tempdir().unwrap();
        let cf_seq_path = dir.path().join("order.cf_seq");
        let json_path = dir.path().join("none.json");
        // File order of segments: 30, 10, 20, 40 (40 appears in no tiling).
        std::fs::write(
            &cf_seq_path,
            "Reference:1_Sequence:a\t20+ N4 10-\nReference:1_Sequence:b\t10+ 30+\n",
        )
        .unwrap();
        let ids = [30u64, 10, 20, 40];
        let mut id_to_rank = FixedMap::with_hasher(crate::hash::fixed_state());
        for (rank, &id) in ids.iter().enumerate() {
            id_to_rank.insert(id, rank as u32);
        }
        let mut tilings = parse_cf_seq(
            cf_seq_path.to_str().unwrap(),
            json_path.to_str().unwrap(),
            5,
            &id_to_rank,
            &[10, 11, 12, 13],
        )
        .unwrap();
        assert_eq!(tilings.counts, [1, 2, 1, 0]);

        // First appearance: 20 (rank 2), 10 (1), 30 (0); then unseen 40 (3).
        let order = first_appearance_order(&tilings, &ids);
        assert_eq!(order, [2, 1, 0, 3]);
        tilings.renumber(&order);
        assert_eq!(tilings.counts, [1, 2, 1, 0]);
        let tiles: Vec<u64> = tilings
            .chunks
            .iter()
            .flat_map(|c| c.tiles.clone())
            .collect();
        // `rank << 1 | is_fw`: a = 20+ N4 10-, b = 10+ 30+.
        assert_eq!(tiles, [1, GAP_TILE | 4, 2, 3, 5]);
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
            seed: 1,
            mphf_build: Default::default(),
            emit_tiny: Some(false),
            tmp_dir: None,
            ram_limit_gib: None,
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
