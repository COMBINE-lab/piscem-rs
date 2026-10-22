//! The full piscem reference index.
//!
//! Assembles the SSHash dictionary, contig table, reference metadata, and
//! optional components (equivalence class table, poison table) into a single
//! loadable/saveable index.
//!
//! Corresponds to the C++ `reference_index` class.

use anyhow::{Context, Result};
use sshash_lib::{Dictionary, Kmer, KmerBits, KmerDictionary};
use std::path::Path;
use tiny_dict::TinyDictionary;
use tracing::{info, warn};

use super::contig_table::{ContigTable, ContigTableLike, EntryEncoding, TinyContigTable};
use super::eq_classes::EqClassMap;
use super::paths::{
    CTAB_EXT, DICT_EXT, ECTAB_EXT, IndexLayout, IndexPaths, POISON_EXT, POISON_JSON_EXT,
    REFINFO_EXT, SIGS_EXT, TCT_EXT, TDCT_EXT,
};
use super::poison_table::PoisonTable;
use super::refinfo::RefInfo;

use crate::mapping::projected_hits::ProjectedHits;
use sshash_lib::LookupResult;

// ---------------------------------------------------------------------------
// RefSigInfo
// ---------------------------------------------------------------------------

/// Reference signature hashes loaded from `{index_prefix}.sigs.json`.
///
/// Matches the C++ `ref_sig_info_t` struct. All four fields are hex-encoded
/// hash strings. The file is optional — if absent, `load()` returns `None`.
#[derive(Debug, Clone)]
pub struct RefSigInfo {
    pub sha256_names: String,
    pub sha256_seqs: String,
    pub sha512_names: String,
    pub sha512_seqs: String,
}

impl RefSigInfo {
    /// From the parsed signature file at `path` (`{prefix}.sigs.json`, read
    /// by [`read_sigs_json`]); `None` if there is none.
    fn from_json(path: &Path, v: Option<&serde_json::Value>) -> Option<Self> {
        let v = v?;
        let get = |key: &str| v.get(key).and_then(|s| s.as_str()).unwrap_or("").to_owned();
        info!("Loaded reference signatures from {}", path.display());
        Some(Self {
            sha256_names: get("sha256_names"),
            sha256_seqs: get("sha256_seqs"),
            sha512_names: get("sha512_names"),
            sha512_seqs: get("sha512_seqs"),
        })
    }
}

/// The parsed signature file, or `None` (with a warning if it exists but is
/// unreadable).
fn read_sigs_json(path: &Path) -> Option<serde_json::Value> {
    if !path.exists() {
        return None;
    }
    let file = match std::fs::File::open(path) {
        Ok(f) => f,
        Err(e) => {
            warn!(
                "Signature file {} exists but could not be opened: {}",
                path.display(),
                e
            );
            return None;
        }
    };
    match serde_json::from_reader(std::io::BufReader::new(file)) {
        Ok(v) => Some(v),
        Err(e) => {
            warn!(
                "Signature file {} could not be parsed: {}",
                path.display(),
                e
            );
            None
        }
    }
}

/// For a [`IndexLayout::Legacy`] index, warn when the core files hold a
/// different number of references than the signature file (`sigs`, as
/// parsed by [`read_sigs_json`]) records: the replaced-suffix paths are
/// shared by every prefix with the same stem, so another index may have
/// overwritten them.
fn check_legacy_ref_count(paths: &IndexPaths, num_refs: usize, sigs: Option<&serde_json::Value>) {
    if paths.layout() != IndexLayout::Legacy {
        return;
    }
    let Some(recorded) = sigs.and_then(|v| v.get("num_records").and_then(|n| n.as_u64())) else {
        return;
    };
    if recorded != num_refs as u64 {
        warn!(
            "{} holds {} references, but {} records {}: the index files were probably \
             overwritten by another index built at a prefix with the same stem, and mapping \
             results will be wrong. Rebuild the index",
            paths.file(REFINFO_EXT).display(),
            num_refs,
            paths.file(SIGS_EXT).display(),
            recorded
        );
    }
}

/// The error for the first index file a load cannot read: `msg`, plus a
/// note when `prefix` ends in a path separator (see
/// [`super::paths::directory_prefix_hint`]).
fn first_file_error(prefix: &Path, msg: String) -> anyhow::Error {
    match super::paths::directory_prefix_hint(prefix) {
        Some(hint) => anyhow::anyhow!("{msg} ({hint})"),
        None => anyhow::anyhow!(msg),
    }
}

/// Find the index at `prefix` (once per load) and report how it was found.
fn resolve_for_load(prefix: &Path) -> Result<IndexPaths> {
    if prefix.file_name().is_none() {
        anyhow::bail!(
            "index prefix {:?} names a directory, not an index: pass the prefix the index \
             was built with (e.g. {:?})",
            prefix,
            prefix.join("index")
        );
    }
    let paths = IndexPaths::resolve(prefix);
    paths.log_resolution();
    Ok(paths)
}

/// Whether both `.tdct` and `.tct` files exist for the index at `prefix`
/// (in whichever layout [`IndexPaths::resolve`] finds it).
pub fn tiny_artifacts_exist(prefix: &Path) -> bool {
    IndexPaths::resolve(prefix).tiny_artifacts_exist()
}

// ---------------------------------------------------------------------------
// ReferenceIndex
// ---------------------------------------------------------------------------

/// The complete piscem reference index.
///
/// Contains:
/// - An SSHash dictionary for k-mer → unitig lookup
/// - A contig table mapping unitigs → packed reference occurrences
/// - Reference metadata (names and lengths)
/// - Derived encoding parameters for decoding contig table entries
///
/// # Loading
///
/// Given a prefix like `/path/to/index`, the following files are expected
/// (each is the full prefix plus an extension, see
/// [`crate::index::paths::index_artifact_path`]):
/// - `/path/to/index.ssi` + `/path/to/index.ssi.mphf` — SSHash dictionary
/// - `/path/to/index.ctab` — contig table
/// - `/path/to/index.refinfo` — reference names and lengths
/// - `/path/to/index.ectab` — equivalence class table (optional)
/// - `/path/to/index.poison` — poison k-mer table (optional)
/// - `/path/to/index.sigs.json` — reference signatures (optional)
///
/// Indexes that piscem-rs 0.11.0 and earlier wrote at a dotted prefix are
/// still loaded, with a warning; see [`IndexPaths::resolve`].
pub struct ReferenceIndex<D = Dictionary, C = ContigTable>
where
    D: KmerDictionary,
    C: ContigTableLike,
{
    /// SSHash compressed k-mer dictionary.
    dict: D,
    /// Contig (unitig) → packed reference occurrence table.
    contig_table: C,
    /// Reference sequence names and lengths.
    ref_info: RefInfo,
    /// Pre-computed entry encoding parameters (derived from contig table).
    encoding: EntryEncoding,
    /// Optional equivalence class table (tile → EC → label entries).
    ec_table: Option<EqClassMap>,
    /// Optional poison k-mer table for mapping specificity.
    poison_table: Option<PoisonTable>,
    /// Optional reference signature hashes from `{prefix}.sigs.json`.
    sig_info: Option<RefSigInfo>,
}

// ---------------------------------------------------------------------------
// File I/O and construction — only for the default (sshash) parameters.
// TinyDictionary-based variants will add their own load/save paths.
// ---------------------------------------------------------------------------
impl ReferenceIndex<Dictionary, ContigTable> {
    /// Load a reference index from files with the given path prefix.
    ///
    /// Loads the dictionary, contig table, and reference info. Optionally
    /// loads the equivalence class table if `load_ec` is true and the
    /// poison table if `load_poison` is true (and their files exist).
    pub fn load(prefix: &Path, load_ec: bool, load_poison: bool) -> Result<Self> {
        let paths = resolve_for_load(prefix)?;

        // 1. Load SSHash dictionary
        let dict_path = paths.file(DICT_EXT);
        info!("Loading SSHash dictionary from {}", dict_path.display());
        let dict = Dictionary::load(&dict_path).map_err(|e| {
            first_file_error(
                prefix,
                format!("failed to load dictionary {}: {e}", dict_path.display()),
            )
        })?;
        info!("  k={}, {} strings", dict.k(), dict.num_strings());

        // 2. Load contig table
        let ctab_path = paths.file(CTAB_EXT);
        info!("Loading contig table from {}", ctab_path.display());
        let ctab_file = std::fs::File::open(&ctab_path)
            .with_context(|| format!("failed to open {}", ctab_path.display()))?;
        let mut ctab_reader = std::io::BufReader::new(ctab_file);
        let contig_table = ContigTable::load(&mut ctab_reader)
            .with_context(|| format!("failed to load {}", ctab_path.display()))?;
        info!(
            "  {} contigs, {} entries, entry_width={} bits",
            contig_table.num_contigs(),
            contig_table.num_entries(),
            contig_table.ref_len_bits() + contig_table.num_ref_bits() + 1,
        );

        // 3. Load reference info
        let refinfo_path = paths.file(REFINFO_EXT);
        info!("Loading reference info from {}", refinfo_path.display());
        let refinfo_file = std::fs::File::open(&refinfo_path)
            .with_context(|| format!("failed to open {}", refinfo_path.display()))?;
        let mut refinfo_reader = std::io::BufReader::new(refinfo_file);
        let ref_info = RefInfo::load(&mut refinfo_reader)
            .with_context(|| format!("failed to load {}", refinfo_path.display()))?;
        info!(
            "  {} references, max_ref_len={}",
            ref_info.num_refs(),
            ref_info.max_ref_len(),
        );
        let sigs_path = paths.file(SIGS_EXT);
        let sigs = read_sigs_json(&sigs_path);
        check_legacy_ref_count(&paths, ref_info.num_refs(), sigs.as_ref());

        // 4. Optionally load equivalence class table
        let ec_table = if load_ec {
            let ectab_path = paths.file(ECTAB_EXT);
            if ectab_path.exists() {
                info!(
                    "Loading equivalence class table from {}",
                    ectab_path.display()
                );
                let ectab_file = std::fs::File::open(&ectab_path)
                    .with_context(|| format!("failed to open {}", ectab_path.display()))?;
                let mut ectab_reader = std::io::BufReader::new(ectab_file);
                let ec = EqClassMap::load(&mut ectab_reader)
                    .with_context(|| format!("failed to load {}", ectab_path.display()))?;
                info!(
                    "  {} tiles, {} ECs, {} label entries",
                    ec.num_tiles(),
                    ec.num_ecs(),
                    ec.num_label_entries(),
                );
                Some(ec)
            } else {
                info!(
                    "No equivalence class table found at {}",
                    ectab_path.display()
                );
                None
            }
        } else {
            None
        };

        // 5. Derive encoding parameters from contig table
        let encoding = contig_table.encoding();
        info!(
            "  ref_shift={}, pos_mask=0x{:x}",
            encoding.ref_shift, encoding.pos_mask,
        );

        // 6. Optionally load poison table
        let poison_table = if load_poison {
            let poison_path = paths.file(POISON_EXT);
            if poison_path.exists() {
                info!("Loading poison table from {}", poison_path.display());
                let poison_file = std::fs::File::open(&poison_path)
                    .with_context(|| format!("failed to open {}", poison_path.display()))?;
                let mut poison_reader = std::io::BufReader::new(poison_file);
                let pt = PoisonTable::load(&mut poison_reader)
                    .with_context(|| format!("failed to load {}", poison_path.display()))?;
                info!(
                    "  {} poison k-mers, {} occurrences",
                    pt.num_poison_kmers(),
                    pt.num_poison_occs(),
                );
                Some(pt)
            } else {
                info!("No poison table found at {}", poison_path.display());
                None
            }
        } else {
            None
        };

        // 7. Optionally load reference signature info
        let sig_info = RefSigInfo::from_json(&sigs_path, sigs.as_ref());

        Ok(Self {
            dict,
            contig_table,
            ref_info,
            encoding,
            ec_table,
            poison_table,
            sig_info,
        })
    }

    /// Save the index to files with the given path prefix.
    pub fn save(&self, prefix: &Path) -> Result<()> {
        let paths = IndexPaths::new(prefix);
        // 1. Save dictionary
        let dict_path = paths.file(DICT_EXT);
        info!("Saving SSHash dictionary to {}", dict_path.display());
        self.dict
            .save(&dict_path)
            .map_err(|e| anyhow::anyhow!("failed to save dictionary: {e}"))?;

        // 2. Save contig table
        let ctab_path = paths.file(CTAB_EXT);
        info!("Saving contig table to {}", ctab_path.display());
        let ctab_file = std::fs::File::create(&ctab_path)
            .with_context(|| format!("failed to create {}", ctab_path.display()))?;
        let mut ctab_writer = std::io::BufWriter::new(ctab_file);
        self.contig_table
            .save(&mut ctab_writer)
            .with_context(|| format!("failed to save {}", ctab_path.display()))?;

        // 3. Save reference info
        let refinfo_path = paths.file(REFINFO_EXT);
        info!("Saving reference info to {}", refinfo_path.display());
        let refinfo_file = std::fs::File::create(&refinfo_path)
            .with_context(|| format!("failed to create {}", refinfo_path.display()))?;
        let mut refinfo_writer = std::io::BufWriter::new(refinfo_file);
        self.ref_info
            .save(&mut refinfo_writer)
            .with_context(|| format!("failed to save {}", refinfo_path.display()))?;

        // 4. Save EC table if present
        if let Some(ref ec) = self.ec_table {
            let ectab_path = paths.file(ECTAB_EXT);
            info!("Saving equivalence class table to {}", ectab_path.display());
            let ectab_file = std::fs::File::create(&ectab_path)
                .with_context(|| format!("failed to create {}", ectab_path.display()))?;
            let mut ectab_writer = std::io::BufWriter::new(ectab_file);
            ec.save(&mut ectab_writer)
                .with_context(|| format!("failed to save {}", ectab_path.display()))?;
        }

        // 5. Save poison table if present
        if let Some(ref pt) = self.poison_table {
            let poison_path = paths.file(POISON_EXT);
            info!("Saving poison table to {}", poison_path.display());
            let poison_file = std::fs::File::create(&poison_path)
                .with_context(|| format!("failed to create {}", poison_path.display()))?;
            let mut poison_writer = std::io::BufWriter::new(poison_file);
            pt.save(&mut poison_writer)
                .with_context(|| format!("failed to save {}", poison_path.display()))?;

            let json_path = paths.file(POISON_JSON_EXT);
            info!("Saving poison stats to {}", json_path.display());
            let json_file = std::fs::File::create(&json_path)
                .with_context(|| format!("failed to create {}", json_path.display()))?;
            let mut json_writer = std::io::BufWriter::new(json_file);
            pt.save_stats_json(&mut json_writer)
                .with_context(|| format!("failed to save {}", json_path.display()))?;
        }

        info!("Index saved to {}", prefix.display());
        Ok(())
    }

    // -----------------------------------------------------------------------
    // Construction (for building a new index)
    // -----------------------------------------------------------------------

    /// Assemble a `ReferenceIndex` from pre-built components.
    pub fn from_parts(
        dict: Dictionary,
        contig_table: ContigTable,
        ref_info: RefInfo,
        ec_table: Option<EqClassMap>,
        poison_table: Option<PoisonTable>,
    ) -> Self {
        let encoding = contig_table.encoding();
        Self {
            dict,
            contig_table,
            ref_info,
            encoding,
            ec_table,
            poison_table,
            sig_info: None,
        }
    }

    /// Convert this sshash-backed index into a TinyDictionary-backed one.
    ///
    /// Builds a `TinyDictionary` from the loaded sshash dictionary and keeps
    /// every other component (contig table, refinfo, EC/poison tables) as-is.
    ///
    /// TODO: eventually persist the TinyDictionary directly (its own `.tdct`
    /// file) so we can skip the sshash → tiny conversion at load time.
    pub fn into_tiny<const K: usize>(self) -> ReferenceIndex<TinyDictionary, ContigTable>
    where
        Kmer<K>: KmerBits,
    {
        let tiny = TinyDictionary::from_sshash::<K>(&self.dict);
        ReferenceIndex {
            dict: tiny,
            contig_table: self.contig_table,
            ref_info: self.ref_info,
            encoding: self.encoding,
            ec_table: self.ec_table,
            poison_table: self.poison_table,
            sig_info: self.sig_info,
        }
    }

    /// Serialize Tiny artifacts (`.tdct` + `.tct`) alongside the existing
    /// sshash `.ssi` / `.ctab` files for this index prefix.
    ///
    /// Converts the loaded sshash dictionary and contig table into their Tiny
    /// counterparts and writes them under `{prefix}.tdct` and `{prefix}.tct`.
    /// The sshash artifacts are untouched; load-time code selects the backend
    /// based on the `--dict` flag.
    pub fn save_tiny_artifacts<const K: usize>(&self, prefix: &Path) -> Result<()>
    where
        Kmer<K>: KmerBits,
    {
        let paths = IndexPaths::new(prefix);
        let tdct_path = paths.file(TDCT_EXT);
        info!(
            "Building and saving TinyDictionary to {}",
            tdct_path.display()
        );
        let tiny_dict = TinyDictionary::from_sshash::<K>(&self.dict);
        tiny_dict
            .save(&tdct_path)
            .map_err(|e| anyhow::anyhow!("failed to save {}: {e}", tdct_path.display()))?;

        let tct_path = paths.file(TCT_EXT);
        info!(
            "Building and saving TinyContigTable to {}",
            tct_path.display()
        );
        let tiny_ctab = TinyContigTable::from_contig_table(&self.contig_table);
        let tct_file = std::fs::File::create(&tct_path)
            .with_context(|| format!("failed to create {}", tct_path.display()))?;
        let mut tct_writer = std::io::BufWriter::new(tct_file);
        tiny_ctab
            .save(&mut tct_writer)
            .with_context(|| format!("failed to save {}", tct_path.display()))?;
        Ok(())
    }

    /// Convert to a fully Tiny-backed index: TinyDictionary + TinyContigTable.
    ///
    /// Builds a `TinyDictionary` from the loaded sshash dictionary and a
    /// `TinyContigTable` (inline-single-ref fast path) from the loaded
    /// `ContigTable`. Other components are preserved.
    pub fn into_tiny_full<const K: usize>(self) -> ReferenceIndex<TinyDictionary, TinyContigTable>
    where
        Kmer<K>: KmerBits,
    {
        let tiny_dict = TinyDictionary::from_sshash::<K>(&self.dict);
        let tiny_ctab = TinyContigTable::from_contig_table(&self.contig_table);
        ReferenceIndex {
            dict: tiny_dict,
            contig_table: tiny_ctab,
            ref_info: self.ref_info,
            encoding: self.encoding,
            ec_table: self.ec_table,
            poison_table: self.poison_table,
            sig_info: self.sig_info,
        }
    }
}

// ---------------------------------------------------------------------------
// Tiny-backed load path: reads `.tdct` + `.tct` directly, sharing refinfo /
// ectab / poison / sigs with the sshash layout.
// ---------------------------------------------------------------------------
impl ReferenceIndex<TinyDictionary, TinyContigTable> {
    /// Load a Tiny-backed index directly from `.tdct` + `.tct` files.
    ///
    /// Expects the same companion files as the sshash loader
    /// (`.refinfo`, optional `.ectab`, `.poison`, `.sigs.json`).
    pub fn load_tiny(prefix: &Path, load_ec: bool, load_poison: bool) -> Result<Self> {
        let paths = resolve_for_load(prefix)?;

        // 1. TinyDictionary
        let tdct_path = paths.file(TDCT_EXT);
        info!("Loading TinyDictionary from {}", tdct_path.display());
        let dict = TinyDictionary::load(&tdct_path).map_err(|e| {
            first_file_error(
                prefix,
                format!("failed to load {}: {e}", tdct_path.display()),
            )
        })?;
        info!("  k={}, {} strings", dict.k(), dict.num_strings());

        // 2. TinyContigTable
        let tct_path = paths.file(TCT_EXT);
        info!("Loading TinyContigTable from {}", tct_path.display());
        let tct_file = std::fs::File::open(&tct_path)
            .with_context(|| format!("failed to open {}", tct_path.display()))?;
        let mut tct_reader = std::io::BufReader::new(tct_file);
        let contig_table = TinyContigTable::load(&mut tct_reader)
            .with_context(|| format!("failed to load {}", tct_path.display()))?;
        info!(
            "  {} contigs, inline-single-ref fraction = {:.2}%",
            contig_table.num_contigs(),
            contig_table.inline_fraction() * 100.0,
        );

        // 3. Reference info
        let refinfo_path = paths.file(REFINFO_EXT);
        let refinfo_file = std::fs::File::open(&refinfo_path)
            .with_context(|| format!("failed to open {}", refinfo_path.display()))?;
        let mut refinfo_reader = std::io::BufReader::new(refinfo_file);
        let ref_info = RefInfo::load(&mut refinfo_reader)
            .with_context(|| format!("failed to load {}", refinfo_path.display()))?;
        let sigs_path = paths.file(SIGS_EXT);
        let sigs = read_sigs_json(&sigs_path);
        check_legacy_ref_count(&paths, ref_info.num_refs(), sigs.as_ref());

        // 4. EC table
        let ec_table = if load_ec {
            let ectab_path = paths.file(ECTAB_EXT);
            if ectab_path.exists() {
                let ectab_file = std::fs::File::open(&ectab_path)
                    .with_context(|| format!("failed to open {}", ectab_path.display()))?;
                let mut ectab_reader = std::io::BufReader::new(ectab_file);
                Some(
                    EqClassMap::load(&mut ectab_reader)
                        .with_context(|| format!("failed to load {}", ectab_path.display()))?,
                )
            } else {
                None
            }
        } else {
            None
        };

        // 5. Poison table
        let poison_table = if load_poison {
            let poison_path = paths.file(POISON_EXT);
            if poison_path.exists() {
                let poison_file = std::fs::File::open(&poison_path)
                    .with_context(|| format!("failed to open {}", poison_path.display()))?;
                let mut poison_reader = std::io::BufReader::new(poison_file);
                Some(
                    PoisonTable::load(&mut poison_reader)
                        .with_context(|| format!("failed to load {}", poison_path.display()))?,
                )
            } else {
                None
            }
        } else {
            None
        };

        let sig_info = RefSigInfo::from_json(&sigs_path, sigs.as_ref());
        let encoding = contig_table.encoding();

        Ok(Self {
            dict,
            contig_table,
            ref_info,
            encoding,
            ec_table,
            poison_table,
            sig_info,
        })
    }
}

// ---------------------------------------------------------------------------
// Generic accessors and query/resolution — work for any KmerDictionary +
// ContigTableLike pair. Today only the default (Dictionary, ContigTable) is
// instantiated; the generics are here so future TinyDictionary / TinyContigTable
// variants can drop in without touching call sites.
// ---------------------------------------------------------------------------
impl<D, C> ReferenceIndex<D, C>
where
    D: KmerDictionary,
    C: ContigTableLike,
{
    // -----------------------------------------------------------------------
    // Accessors
    // -----------------------------------------------------------------------

    /// K-mer size used by the dictionary.
    #[inline]
    pub fn k(&self) -> usize {
        self.dict.k()
    }

    /// Minimizer length used by the dictionary.
    #[inline]
    pub fn m(&self) -> usize {
        self.dict.m()
    }

    /// The SSHash dictionary.
    #[inline]
    pub fn dict(&self) -> &D {
        &self.dict
    }

    /// The contig (unitig) occurrence table.
    #[inline]
    pub fn contig_table(&self) -> &C {
        &self.contig_table
    }

    /// Reference sequence metadata.
    #[inline]
    pub fn ref_info(&self) -> &RefInfo {
        &self.ref_info
    }

    /// Number of reference sequences.
    #[inline]
    pub fn num_refs(&self) -> usize {
        self.ref_info.num_refs()
    }

    /// Name of reference `i`.
    #[inline]
    pub fn ref_name(&self, i: usize) -> &str {
        self.ref_info.name(i)
    }

    /// Length of reference `i`.
    #[inline]
    pub fn ref_len(&self, i: usize) -> u64 {
        self.ref_info.len(i)
    }

    /// Pre-computed entry encoding parameters for decoding contig table entries.
    #[inline]
    pub fn encoding(&self) -> EntryEncoding {
        self.encoding
    }

    /// Whether an equivalence class table is available.
    #[inline]
    pub fn has_ec_table(&self) -> bool {
        self.ec_table.is_some()
    }

    /// The equivalence class table, if loaded.
    #[inline]
    pub fn ec_table(&self) -> Option<&EqClassMap> {
        self.ec_table.as_ref()
    }

    /// Whether a poison k-mer table is available.
    #[inline]
    pub fn has_poison_table(&self) -> bool {
        self.poison_table.is_some()
    }

    /// The poison k-mer table, if loaded.
    #[inline]
    pub fn poison_table(&self) -> Option<&PoisonTable> {
        self.poison_table.as_ref()
    }

    /// Reference signature hashes, if a `.sigs.json` file was found alongside the index.
    #[inline]
    pub fn ref_sig_info(&self) -> Option<&RefSigInfo> {
        self.sig_info.as_ref()
    }

    /// Number of unitigs (contigs) in the index.
    #[inline]
    pub fn num_contigs(&self) -> usize {
        self.contig_table.num_contigs()
    }

    // -----------------------------------------------------------------------
    // Query / resolution
    // -----------------------------------------------------------------------

    /// Resolve an sshash-rs `LookupResult` against the contig table.
    ///
    /// Returns `None` if the k-mer was not found in the dictionary.
    /// Returns `Some(hits)` where `hits.is_empty()` is possible if the k-mer
    /// was found but the contig has no reference occurrences.
    pub fn resolve_lookup(&self, result: &LookupResult) -> Option<ProjectedHits<'_>> {
        if !result.is_found() {
            return None;
        }
        let k = self.k() as u32;
        let contig_id = result.string_id as u32;
        let contig_pos = result.kmer_id_in_string as u32;
        let contig_len = result.string_length() as u32;
        let is_forward = result.kmer_orientation == 1;
        let global_pos = result.kmer_id + result.string_id * (k as u64 - 1);
        let ref_range = self.contig_table.contig_entries(result.string_id);
        Some(ProjectedHits::new(
            contig_id, contig_pos, is_forward, contig_len, global_pos, k, ref_range,
        ))
    }
}

impl<D, C> std::fmt::Debug for ReferenceIndex<D, C>
where
    D: KmerDictionary,
    C: ContigTableLike,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ReferenceIndex")
            .field("k", &self.k())
            .field("num_contigs", &self.num_contigs())
            .field("num_refs", &self.num_refs())
            .field("has_ec_table", &self.has_ec_table())
            .field("has_poison_table", &self.has_poison_table())
            .field("encoding", &self.encoding)
            .finish()
    }
}
