//! The full piscem reference index.
//!
//! Assembles the SSHash dictionary, contig table, reference metadata, and
//! optional components (equivalence class table, poison table) into a single
//! loadable/saveable index.
//!
//! Corresponds to the C++ `reference_index` class.

use anyhow::{Context, Result};
use sshash_lib::Dictionary;
use std::path::{Path, PathBuf};
use tracing::info;

use super::contig_table::{ContigTable, EntryEncoding};
use super::eq_classes::EqClassMap;
use super::poison_table::PoisonTable;
use super::refinfo::RefInfo;

// ---------------------------------------------------------------------------
// File extension helpers
// ---------------------------------------------------------------------------

/// File extension for the sshash dictionary (base path for `.ssi` + `.ssi.mphf`).
const DICT_EXT: &str = "ssi";
/// File extension for the contig table.
const CTAB_EXT: &str = "ctab";
/// File extension for reference info.
const REFINFO_EXT: &str = "refinfo";
/// File extension for the equivalence class table (optional).
const ECTAB_EXT: &str = "ectab";
/// File extension for the poison table (optional).
const POISON_EXT: &str = "poison";
/// File extension for the poison stats JSON (optional).
const POISON_JSON_EXT: &str = "poison.json";

fn with_ext(prefix: &Path, ext: &str) -> PathBuf {
    let mut p = prefix.to_path_buf();
    p.set_extension(ext);
    p
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
/// Given a prefix like `/path/to/index`, the following files are expected:
/// - `/path/to/index.ssi` + `/path/to/index.ssi.mphf` — SSHash dictionary
/// - `/path/to/index.ctab` — contig table
/// - `/path/to/index.refinfo` — reference names and lengths
/// - `/path/to/index.ectab` — equivalence class table (optional)
/// - `/path/to/index.poison` — poison k-mer table (optional)
pub struct ReferenceIndex {
    /// SSHash compressed k-mer dictionary.
    dict: Dictionary,
    /// Contig (unitig) → packed reference occurrence table.
    contig_table: ContigTable,
    /// Reference sequence names and lengths.
    ref_info: RefInfo,
    /// Pre-computed entry encoding parameters (derived from contig table).
    encoding: EntryEncoding,
    /// Optional equivalence class table (tile → EC → label entries).
    ec_table: Option<EqClassMap>,
    /// Optional poison k-mer table for mapping specificity.
    poison_table: Option<PoisonTable>,
}

impl ReferenceIndex {
    /// Load a reference index from files with the given path prefix.
    ///
    /// Loads the dictionary, contig table, and reference info. Optionally
    /// loads the equivalence class table if `load_ec` is true and the
    /// poison table if `load_poison` is true (and their files exist).
    pub fn load(prefix: &Path, load_ec: bool, load_poison: bool) -> Result<Self> {
        // 1. Load SSHash dictionary
        let dict_path = with_ext(prefix, DICT_EXT);
        info!("Loading SSHash dictionary from {}", dict_path.display());
        let dict = Dictionary::load(&dict_path)
            .map_err(|e| anyhow::anyhow!("failed to load dictionary: {e}"))?;
        info!(
            "  k={}, {} strings, canonical={}",
            dict.k(),
            dict.num_strings(),
            dict.canonical()
        );

        // 2. Load contig table
        let ctab_path = with_ext(prefix, CTAB_EXT);
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
        let refinfo_path = with_ext(prefix, REFINFO_EXT);
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

        // 4. Optionally load equivalence class table
        let ec_table = if load_ec {
            let ectab_path = with_ext(prefix, ECTAB_EXT);
            if ectab_path.exists() {
                info!("Loading equivalence class table from {}", ectab_path.display());
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
                info!("No equivalence class table found at {}", ectab_path.display());
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
            let poison_path = with_ext(prefix, POISON_EXT);
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

        Ok(Self {
            dict,
            contig_table,
            ref_info,
            encoding,
            ec_table,
            poison_table,
        })
    }

    /// Save the index to files with the given path prefix.
    pub fn save(&self, prefix: &Path) -> Result<()> {
        // 1. Save dictionary
        let dict_path = with_ext(prefix, DICT_EXT);
        info!("Saving SSHash dictionary to {}", dict_path.display());
        self.dict
            .save(&dict_path)
            .map_err(|e| anyhow::anyhow!("failed to save dictionary: {e}"))?;

        // 2. Save contig table
        let ctab_path = with_ext(prefix, CTAB_EXT);
        info!("Saving contig table to {}", ctab_path.display());
        let ctab_file = std::fs::File::create(&ctab_path)
            .with_context(|| format!("failed to create {}", ctab_path.display()))?;
        let mut ctab_writer = std::io::BufWriter::new(ctab_file);
        self.contig_table
            .save(&mut ctab_writer)
            .with_context(|| format!("failed to save {}", ctab_path.display()))?;

        // 3. Save reference info
        let refinfo_path = with_ext(prefix, REFINFO_EXT);
        info!("Saving reference info to {}", refinfo_path.display());
        let refinfo_file = std::fs::File::create(&refinfo_path)
            .with_context(|| format!("failed to create {}", refinfo_path.display()))?;
        let mut refinfo_writer = std::io::BufWriter::new(refinfo_file);
        self.ref_info
            .save(&mut refinfo_writer)
            .with_context(|| format!("failed to save {}", refinfo_path.display()))?;

        // 4. Save EC table if present
        if let Some(ref ec) = self.ec_table {
            let ectab_path = with_ext(prefix, ECTAB_EXT);
            info!("Saving equivalence class table to {}", ectab_path.display());
            let ectab_file = std::fs::File::create(&ectab_path)
                .with_context(|| format!("failed to create {}", ectab_path.display()))?;
            let mut ectab_writer = std::io::BufWriter::new(ectab_file);
            ec.save(&mut ectab_writer)
                .with_context(|| format!("failed to save {}", ectab_path.display()))?;
        }

        // 5. Save poison table if present
        if let Some(ref pt) = self.poison_table {
            let poison_path = with_ext(prefix, POISON_EXT);
            info!("Saving poison table to {}", poison_path.display());
            let poison_file = std::fs::File::create(&poison_path)
                .with_context(|| format!("failed to create {}", poison_path.display()))?;
            let mut poison_writer = std::io::BufWriter::new(poison_file);
            pt.save(&mut poison_writer)
                .with_context(|| format!("failed to save {}", poison_path.display()))?;

            let json_path = with_ext(prefix, POISON_JSON_EXT);
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
    // Accessors
    // -----------------------------------------------------------------------

    /// K-mer size used by the dictionary.
    #[inline]
    pub fn k(&self) -> usize {
        self.dict.k()
    }

    /// The SSHash dictionary.
    #[inline]
    pub fn dict(&self) -> &Dictionary {
        &self.dict
    }

    /// The contig (unitig) occurrence table.
    #[inline]
    pub fn contig_table(&self) -> &ContigTable {
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

    /// Number of unitigs (contigs) in the index.
    #[inline]
    pub fn num_contigs(&self) -> usize {
        self.contig_table.num_contigs()
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
        }
    }
}

impl std::fmt::Debug for ReferenceIndex {
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
