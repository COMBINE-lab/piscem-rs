//! Poison k-mer table — improves mapping specificity by identifying k-mers
//! from decoy sequences that sit at the boundary between present/absent on
//! reference unitigs.
//!
//! During mapping, if a read's hits contain a poison k-mer, the mapping is
//! discarded. This avoids false-positive pseudoalignments without requiring
//! alignment scoring.
//!
//! The table maps canonical k-mer values to their occurrence positions on
//! unitigs. The occurrences are stored in a flat array with an offset vector
//! for efficient range access.
//!
//! Corresponds to the C++ `poison_table` class in `poison_table.hpp`.

use ahash::RandomState;
use anyhow::{Context, Result, bail};
use std::collections::HashMap;
use std::io::{BufReader, BufWriter, Read, Write};
use std::sync::OnceLock;
use tracing::info;

use crate::mapping::kmer_value::CanonicalKmer;

// ---------------------------------------------------------------------------
// Fixed hash state for deterministic AHashMap
// ---------------------------------------------------------------------------

/// Create a deterministic `RandomState` for the poison map.
/// Fixed seeds ensure identical iteration order across runs.
fn fixed_hash_state() -> RandomState {
    RandomState::with_seeds(
        0x517cc1b727220a95,
        0x6c62272e07bb0142,
        0x62b821756295c58d,
        0x30b4d5bd83fac2e9,
    )
}

/// Type alias for the poison k-mer hash map.
type PoisonMap = HashMap<CanonicalKmer, u64, RandomState>;

// ---------------------------------------------------------------------------
// Read-side pre-filter
// ---------------------------------------------------------------------------

/// A blocked Bloom filter over the poison k-mers, keyed the way a read is
/// scanned rather than the way the table is keyed.
///
/// Poison scanning visits every k-mer of every mapped read, and essentially
/// every one of them is absent from the table: poison k-mers are only the
/// decoy k-mers on present/absent edges (38 of them for the 10x Flex human
/// `included=FALSE` probes). Answering each visit with a canonical k-mer
/// (rolling a reverse-complement word, taking the min) and an ahash +
/// SwissTable probe cost more than the rest of mapping put together.
///
/// So the filter is keyed by the *forward* k-mer word in the hit searcher's
/// read encoding — base `i` of the window in bits `2i..2i+2`, `(b >> 1) & 3`
/// per base (`A=0 C=1 T=2 G=3`) — and holds *both* orientations of
/// every poison k-mer. A read scan then rolls one word per base and pays one
/// multiply and one L1-resident load per k-mer. Only filter hits (real poison
/// k-mers plus a small false-positive rate) need the canonical k-mer and the
/// exact table query.
///
/// Each key sets two bits within one `u64` word; the filter is sized at ~16
/// bits per inserted word, and at least 32 Kbit (4 KiB, L1-resident) so small
/// tables get a negligible false-positive rate. It is derived from the keys and
/// `k`, so it is built on first use and is not part of the serialized table.
pub(crate) struct ReadKmerFilter {
    words: Vec<u64>,
    /// `64 - log2(words.len())`: shift that maps a hash to a word index.
    word_shift: u32,
    k: u32,
}

impl ReadKmerFilter {
    const MIN_LOG2_WORDS: u32 = 9;
    const BITS_PER_ENTRY: usize = 16;

    fn build<'a>(keys: impl ExactSizeIterator<Item = &'a CanonicalKmer>, k: u32) -> Self {
        assert!(
            (1..=31).contains(&k),
            "poison k-mers support 1 <= k <= 31, got {k}"
        );
        let n_entries = 2 * keys.len();
        let n_words = (n_entries * Self::BITS_PER_ENTRY / 64)
            .next_power_of_two()
            .max(1 << Self::MIN_LOG2_WORDS);
        let mut f = Self {
            words: vec![0u64; n_words],
            word_shift: 64 - n_words.trailing_zeros(),
            k,
        };
        for &km in keys {
            let (fw, rc) = Self::read_words_of(km, k);
            for word in [fw, rc] {
                let (w, bits) = f.slot(word);
                f.words[w] |= bits;
            }
        }
        f
    }

    /// The two read-encoding words a canonical k-mer can appear as.
    ///
    /// `km` is packed first-base-most-significant with `A=0 C=1 G=2 T=3`.
    fn read_words_of(km: CanonicalKmer, k: u32) -> (u64, u64) {
        // ACGT code -> read code (A=0 C=1 T=2 G=3): swap G and T.
        let to_read = |c: u64| c ^ (c >> 1);
        let v = km.as_u64();
        let (mut fw, mut rc) = (0u64, 0u64);
        for i in 0..k {
            let c = (v >> (2 * (k - 1 - i))) & 3;
            fw |= to_read(c) << (2 * i);
            rc |= to_read(c ^ 3) << (2 * (k - 1 - i));
        }
        (fw, rc)
    }

    /// Word index and two-bit mask for a read-encoded k-mer word.
    #[inline(always)]
    fn slot(&self, word: u64) -> (usize, u64) {
        // Multiplicative hashing: the high product bits depend on all of the
        // (at most 62) input bits, so take the word index from the top and
        // the two bit positions from just below it.
        let h = word.wrapping_mul(0x9E37_79B9_7F4A_7C15);
        let w = (h >> self.word_shift) as usize;
        let b1 = (h >> (self.word_shift - 6)) & 63;
        let b2 = (h >> (self.word_shift - 12)) & 63;
        (w, (1u64 << b1) | (1u64 << b2))
    }

    /// `false` means the k-mer whose read-encoded forward word is `word` is
    /// definitely not a poison k-mer (in either orientation).
    #[inline(always)]
    pub(crate) fn may_contain(&self, word: u64) -> bool {
        let (w, bits) = self.slot(word);
        // SAFETY: `w < words.len()`: `h >> word_shift` keeps only the top
        // log2(words.len()) bits of `h`.
        let filter_word = unsafe { *self.words.get_unchecked(w) };
        filter_word & bits == bits
    }
}

// ---------------------------------------------------------------------------
// Occurrence types
// ---------------------------------------------------------------------------

/// A poison k-mer occurrence on a unitig (stored form, no k-mer value).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PoisonOcc {
    /// Unitig (contig) where this poison k-mer occurs.
    pub unitig_id: u32,
    /// Position within the unitig.
    pub unitig_pos: u32,
}

/// A labeled poison occurrence including the canonical k-mer value.
///
/// Used during construction. The `canonical_kmer` field is dropped when
/// building the final `PoisonTable` (the k-mer becomes the hash map key).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LabeledPoisonOcc {
    /// The canonical form of the poison k-mer.
    pub canonical_kmer: CanonicalKmer,
    /// Unitig (contig) where this poison k-mer occurs.
    pub unitig_id: u32,
    /// Position within the unitig.
    pub unitig_pos: u32,
}

// ---------------------------------------------------------------------------
// PoisonTable
// ---------------------------------------------------------------------------

/// Magic bytes for the piscem-rs poison table file format.
const POISON_MAGIC: &[u8; 8] = b"PPOIS01\0";

/// Poison k-mer table for improving mapping specificity.
///
/// Maps canonical k-mer values to their occurrence positions on reference
/// unitigs. For each k-mer, a contiguous range in the `poison_occs` array
/// stores all (unitig_id, unitig_pos) pairs where the poison k-mer occurs.
///
/// Query workflow:
/// 1. Look up `idx = poison_map[kmer]`
/// 2. Occurrences span `poison_occs[offsets[idx]..offsets[idx+1]]`
pub struct PoisonTable {
    /// Hash map: canonical k-mer → index into `offsets`.
    poison_map: PoisonMap,
    /// Cumulative boundaries into `poison_occs`.
    /// For k-mer with map value `i`, occurrences span
    /// `poison_occs[offsets[i]..offsets[i+1]]`.
    offsets: Vec<u32>,
    /// Flat array of poison k-mer occurrences.
    poison_occs: Vec<PoisonOcc>,
    /// Read-side pre-filter, built on first use; see [`ReadKmerFilter`].
    read_filter: OnceLock<ReadKmerFilter>,
}

impl PoisonTable {
    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------

    /// Build a `PoisonTable` from a collection of labeled poison k-mer
    /// occurrences.
    ///
    /// The input is consumed and deduplicated. Matches the C++
    /// `poison_table::build_from_occs()` algorithm:
    /// 1. Sort by (canonical_kmer, unitig_id, unitig_pos)
    /// 2. Remove duplicates
    /// 3. Build offset + map + flat occurrence structures
    pub fn build_from_occs(mut occs: Vec<LabeledPoisonOcc>) -> Result<Self> {
        // Sort by (canonical_kmer, unitig_id, unitig_pos)
        occs.sort_unstable_by(|a, b| {
            a.canonical_kmer
                .cmp(&b.canonical_kmer)
                .then(a.unitig_id.cmp(&b.unitig_id))
                .then(a.unitig_pos.cmp(&b.unitig_pos))
        });

        // Remove duplicates
        occs.dedup();

        info!("Total number of distinct poison k-mer occs: {}", occs.len());

        if occs.is_empty() {
            return Ok(Self {
                read_filter: OnceLock::new(),
                poison_map: HashMap::with_hasher(fixed_hash_state()),
                offsets: Vec::new(),
                poison_occs: Vec::new(),
            });
        }

        // Check that occurrence count fits in u32 offsets
        if occs.len() > u32::MAX as usize {
            bail!(
                "Too many poison k-mer occurrences ({}) for u32 offset vector (max {})",
                occs.len(),
                u32::MAX
            );
        }

        // Build the poison map and offset vector
        let mut poison_map: PoisonMap =
            HashMap::with_capacity_and_hasher(occs.len(), fixed_hash_state());
        let mut offsets: Vec<u32> = Vec::with_capacity(occs.len() + 1);

        let mut range_start_kmer = occs[0].canonical_kmer;
        poison_map.insert(range_start_kmer, offsets.len() as u64);
        offsets.push(0);

        for (i, occ) in occs.iter().enumerate().skip(1) {
            if occ.canonical_kmer != range_start_kmer {
                // New k-mer range starts at position i
                offsets.push(i as u32);
                poison_map.insert(occ.canonical_kmer, (offsets.len() - 1) as u64);
                range_start_kmer = occ.canonical_kmer;
            }
        }

        // Final sentinel: end of last range
        offsets.push(occs.len() as u32);
        offsets.shrink_to_fit();

        // Copy occurrences, dropping the canonical_kmer field
        let poison_occs: Vec<PoisonOcc> = occs
            .iter()
            .map(|o| PoisonOcc {
                unitig_id: o.unitig_id,
                unitig_pos: o.unitig_pos,
            })
            .collect();

        let max_range = offsets.windows(2).map(|w| w[1] - w[0]).max().unwrap_or(0);
        info!(
            "[poison_table]: The most frequently occurring poison k-mer appeared in {} distinct unitig positions.",
            max_range
        );

        Ok(Self {
            read_filter: OnceLock::new(),
            poison_map,
            offsets,
            poison_occs,
        })
    }

    // -----------------------------------------------------------------------
    // Query methods
    // -----------------------------------------------------------------------

    /// Whether this poison table is empty (contains no poison k-mers).
    #[inline]
    pub fn empty(&self) -> bool {
        self.poison_map.is_empty()
    }

    /// Number of distinct poison k-mers.
    #[inline]
    pub fn num_poison_kmers(&self) -> usize {
        self.poison_map.len()
    }

    /// Total number of poison k-mer occurrences.
    #[inline]
    pub fn num_poison_occs(&self) -> usize {
        self.poison_occs.len()
    }

    /// Returns `true` if `km` is a poison k-mer.
    #[inline]
    pub fn key_exists(&self, km: CanonicalKmer) -> bool {
        self.poison_map.contains_key(&km)
    }

    /// Returns `true` if poison k-mer `km` occurs on unitig `unitig_id`.
    #[inline]
    pub fn key_occurs_in_unitig(&self, km: CanonicalKmer, unitig_id: u32) -> bool {
        if let Some(&idx) = self.poison_map.get(&km) {
            let start = self.offsets[idx as usize] as usize;
            let end = self.offsets[idx as usize + 1] as usize;
            self.poison_occs[start..end]
                .iter()
                .any(|occ| occ.unitig_id == unitig_id)
        } else {
            false
        }
    }

    /// Returns `true` if poison k-mer `km` occurs on unitig `unitig_id`
    /// between positions `lb` and `ub` (inclusive).
    #[inline]
    pub fn key_occurs_in_unitig_between(
        &self,
        km: CanonicalKmer,
        unitig_id: u32,
        lb: u32,
        ub: u32,
    ) -> bool {
        if let Some(&idx) = self.poison_map.get(&km) {
            let start = self.offsets[idx as usize] as usize;
            let end = self.offsets[idx as usize + 1] as usize;
            self.poison_occs[start..end].iter().any(|occ| {
                occ.unitig_id == unitig_id && occ.unitig_pos >= lb && occ.unitig_pos <= ub
            })
        } else {
            false
        }
    }

    /// Iterator over all canonical k-mer keys in the poison table.
    pub fn keys(&self) -> impl Iterator<Item = &CanonicalKmer> {
        self.poison_map.keys()
    }

    /// Maximum number of occurrences for any single poison k-mer.
    pub fn max_poison_occ(&self) -> u32 {
        self.offsets
            .windows(2)
            .map(|w| w[1] - w[0])
            .max()
            .unwrap_or(0)
    }

    /// The read-side pre-filter for k-mer length `k`, built on first use.
    ///
    /// An index has a single `k`, so the filter is built once and shared by
    /// all mapping threads.
    #[inline]
    pub(crate) fn read_filter(&self, k: u32) -> &ReadKmerFilter {
        let f = self
            .read_filter
            .get_or_init(|| ReadKmerFilter::build(self.poison_map.keys(), k));
        debug_assert_eq!(f.k, k, "poison read filter requested for a different k");
        f
    }

    // -----------------------------------------------------------------------
    // Serialization
    // -----------------------------------------------------------------------

    /// Serialize the poison table to a writer.
    ///
    /// Format:
    /// ```text
    /// [magic: 8 bytes "PPOIS01\0"]
    /// [num_kmers: u64 LE]
    /// [num_offsets: u64 LE]
    /// [offsets: num_offsets × u32 LE]
    /// [num_occs: u64 LE]
    /// [poison_occs: num_occs × {unitig_id: u32 LE, unitig_pos: u32 LE}]
    /// [map_entries: num_kmers × {canonical_kmer: u64 LE, offset_idx: u64 LE}]
    /// ```
    pub fn save<W: Write>(&self, writer: &mut W) -> Result<()> {
        let mut w = BufWriter::new(writer);

        w.write_all(POISON_MAGIC)?;

        let num_kmers = self.poison_map.len() as u64;
        let num_offsets = self.offsets.len() as u64;
        let num_occs = self.poison_occs.len() as u64;

        // Header counts
        w.write_all(&num_kmers.to_le_bytes())?;

        // Offsets
        w.write_all(&num_offsets.to_le_bytes())?;
        for &off in &self.offsets {
            w.write_all(&off.to_le_bytes())?;
        }

        // Occurrences
        w.write_all(&num_occs.to_le_bytes())?;
        for occ in &self.poison_occs {
            w.write_all(&occ.unitig_id.to_le_bytes())?;
            w.write_all(&occ.unitig_pos.to_le_bytes())?;
        }

        // Hash map entries (sorted by key for deterministic output)
        let mut entries: Vec<(&CanonicalKmer, &u64)> = self.poison_map.iter().collect();
        entries.sort_unstable_by_key(|&(k, _)| *k);
        for (&kmer, &idx) in entries {
            w.write_all(&kmer.as_u64().to_le_bytes())?;
            w.write_all(&idx.to_le_bytes())?;
        }

        w.flush()?;
        Ok(())
    }

    /// Deserialize a poison table from a reader.
    pub fn load<R: Read>(reader: &mut R) -> Result<Self> {
        let mut r = BufReader::new(reader);

        let mut magic = [0u8; 8];
        r.read_exact(&mut magic)
            .context("failed to read poison table magic")?;
        if magic != *POISON_MAGIC {
            bail!(
                "invalid poison table magic: expected {:?}, got {:?}",
                POISON_MAGIC,
                magic
            );
        }

        let num_kmers = read_u64_le(&mut r).context("failed to read num_kmers")? as usize;

        // Offsets
        let num_offsets = read_u64_le(&mut r).context("failed to read num_offsets")? as usize;
        let mut offsets = Vec::with_capacity(num_offsets);
        for _ in 0..num_offsets {
            offsets.push(read_u32_le(&mut r)?);
        }

        // Occurrences
        let num_occs = read_u64_le(&mut r).context("failed to read num_occs")? as usize;
        let mut poison_occs = Vec::with_capacity(num_occs);
        for _ in 0..num_occs {
            let unitig_id = read_u32_le(&mut r)?;
            let unitig_pos = read_u32_le(&mut r)?;
            poison_occs.push(PoisonOcc {
                unitig_id,
                unitig_pos,
            });
        }

        // Hash map entries
        let mut poison_map: PoisonMap =
            HashMap::with_capacity_and_hasher(num_kmers, fixed_hash_state());
        for _ in 0..num_kmers {
            let kmer = CanonicalKmer::new(read_u64_le(&mut r)?);
            let idx = read_u64_le(&mut r)?;
            poison_map.insert(kmer, idx);
        }

        Ok(Self {
            read_filter: OnceLock::new(),
            poison_map,
            offsets,
            poison_occs,
        })
    }

    /// Save a JSON stats file alongside the poison table.
    pub fn save_stats_json<W: Write>(&self, writer: &mut W) -> Result<()> {
        let stats = serde_json::json!({
            "num_poison_kmers": self.num_poison_kmers(),
            "num_poison_occs": self.num_poison_occs(),
            "max_poison_occ": self.max_poison_occ(),
        });
        serde_json::to_writer_pretty(writer, &stats)
            .context("failed to write poison stats JSON")?;
        Ok(())
    }
}

impl std::fmt::Debug for PoisonTable {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PoisonTable")
            .field("num_poison_kmers", &self.num_poison_kmers())
            .field("num_poison_occs", &self.num_poison_occs())
            .field("max_poison_occ", &self.max_poison_occ())
            .finish()
    }
}

// ---------------------------------------------------------------------------
// I/O helpers
// ---------------------------------------------------------------------------

fn read_u64_le<R: Read>(reader: &mut R) -> std::io::Result<u64> {
    let mut buf = [0u8; 8];
    reader.read_exact(&mut buf)?;
    Ok(u64::from_le_bytes(buf))
}

fn read_u32_le<R: Read>(reader: &mut R) -> std::io::Result<u32> {
    let mut buf = [0u8; 4];
    reader.read_exact(&mut buf)?;
    Ok(u32::from_le_bytes(buf))
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn make_test_occs() -> Vec<LabeledPoisonOcc> {
        vec![
            // K-mer 100: occurs on unitig 0 at pos 5, unitig 1 at pos 10
            LabeledPoisonOcc {
                canonical_kmer: CanonicalKmer::new(100),
                unitig_id: 0,
                unitig_pos: 5,
            },
            LabeledPoisonOcc {
                canonical_kmer: CanonicalKmer::new(100),
                unitig_id: 1,
                unitig_pos: 10,
            },
            // K-mer 200: occurs on unitig 2 at pos 0
            LabeledPoisonOcc {
                canonical_kmer: CanonicalKmer::new(200),
                unitig_id: 2,
                unitig_pos: 0,
            },
            // K-mer 300: occurs on unitig 0 at pos 20, pos 25, pos 30
            LabeledPoisonOcc {
                canonical_kmer: CanonicalKmer::new(300),
                unitig_id: 0,
                unitig_pos: 20,
            },
            LabeledPoisonOcc {
                canonical_kmer: CanonicalKmer::new(300),
                unitig_id: 0,
                unitig_pos: 25,
            },
            LabeledPoisonOcc {
                canonical_kmer: CanonicalKmer::new(300),
                unitig_id: 0,
                unitig_pos: 30,
            },
        ]
    }

    #[test]
    fn test_build_from_occs_basic() {
        let occs = make_test_occs();
        let table = PoisonTable::build_from_occs(occs).unwrap();

        assert_eq!(table.num_poison_kmers(), 3);
        assert_eq!(table.num_poison_occs(), 6);
        assert!(!table.empty());
    }

    #[test]
    fn test_build_from_occs_deduplication() {
        let mut occs = make_test_occs();
        // Add duplicates
        occs.push(LabeledPoisonOcc {
            canonical_kmer: CanonicalKmer::new(100),
            unitig_id: 0,
            unitig_pos: 5,
        });
        occs.push(LabeledPoisonOcc {
            canonical_kmer: CanonicalKmer::new(200),
            unitig_id: 2,
            unitig_pos: 0,
        });

        let table = PoisonTable::build_from_occs(occs).unwrap();

        // Duplicates should be removed
        assert_eq!(table.num_poison_kmers(), 3);
        assert_eq!(table.num_poison_occs(), 6);
    }

    #[test]
    fn test_key_exists() {
        let table = PoisonTable::build_from_occs(make_test_occs()).unwrap();

        assert!(table.key_exists(CanonicalKmer::new(100)));
        assert!(table.key_exists(CanonicalKmer::new(200)));
        assert!(table.key_exists(CanonicalKmer::new(300)));
        assert!(!table.key_exists(CanonicalKmer::new(999)));
        assert!(!table.key_exists(CanonicalKmer::new(0)));
    }

    #[test]
    fn test_key_occurs_in_unitig() {
        let table = PoisonTable::build_from_occs(make_test_occs()).unwrap();

        // K-mer 100 occurs on unitigs 0 and 1
        assert!(table.key_occurs_in_unitig(CanonicalKmer::new(100), 0));
        assert!(table.key_occurs_in_unitig(CanonicalKmer::new(100), 1));
        assert!(!table.key_occurs_in_unitig(CanonicalKmer::new(100), 2));

        // K-mer 200 occurs only on unitig 2
        assert!(!table.key_occurs_in_unitig(CanonicalKmer::new(200), 0));
        assert!(table.key_occurs_in_unitig(CanonicalKmer::new(200), 2));

        // K-mer 300 occurs only on unitig 0
        assert!(table.key_occurs_in_unitig(CanonicalKmer::new(300), 0));
        assert!(!table.key_occurs_in_unitig(CanonicalKmer::new(300), 1));

        // Non-existent k-mer
        assert!(!table.key_occurs_in_unitig(CanonicalKmer::new(999), 0));
    }

    #[test]
    fn test_key_occurs_in_unitig_between() {
        let table = PoisonTable::build_from_occs(make_test_occs()).unwrap();

        // K-mer 300 on unitig 0: positions 20, 25, 30
        let km300 = CanonicalKmer::new(300);
        assert!(table.key_occurs_in_unitig_between(km300, 0, 20, 30));
        assert!(table.key_occurs_in_unitig_between(km300, 0, 20, 20));
        assert!(table.key_occurs_in_unitig_between(km300, 0, 25, 25));
        assert!(table.key_occurs_in_unitig_between(km300, 0, 24, 26));
        assert!(!table.key_occurs_in_unitig_between(km300, 0, 21, 24));
        assert!(!table.key_occurs_in_unitig_between(km300, 0, 26, 29));
        assert!(!table.key_occurs_in_unitig_between(km300, 0, 31, 40));

        // Wrong unitig
        assert!(!table.key_occurs_in_unitig_between(km300, 1, 20, 30));

        // Non-existent k-mer
        assert!(!table.key_occurs_in_unitig_between(CanonicalKmer::new(999), 0, 0, 100));
    }

    #[test]
    fn test_empty_table() {
        let table = PoisonTable::build_from_occs(Vec::new()).unwrap();

        assert!(table.empty());
        assert_eq!(table.num_poison_kmers(), 0);
        assert_eq!(table.num_poison_occs(), 0);
        assert_eq!(table.max_poison_occ(), 0);
        assert!(!table.key_exists(CanonicalKmer::new(42)));
        assert!(!table.key_occurs_in_unitig(CanonicalKmer::new(42), 0));
        assert!(!table.key_occurs_in_unitig_between(CanonicalKmer::new(42), 0, 0, 100));
    }

    #[test]
    fn test_max_poison_occ() {
        let table = PoisonTable::build_from_occs(make_test_occs()).unwrap();

        // K-mer 300 has 3 occurrences (the max)
        assert_eq!(table.max_poison_occ(), 3);
    }

    #[test]
    fn test_serialization_roundtrip() {
        let table = PoisonTable::build_from_occs(make_test_occs()).unwrap();

        // Serialize
        let mut buf = Vec::new();
        table.save(&mut buf).unwrap();

        // Deserialize
        let table2 = PoisonTable::load(&mut &buf[..]).unwrap();

        // Verify structure
        assert_eq!(table2.num_poison_kmers(), table.num_poison_kmers());
        assert_eq!(table2.num_poison_occs(), table.num_poison_occs());
        assert_eq!(table2.max_poison_occ(), table.max_poison_occ());

        // Verify all queries produce the same results
        for &raw in &[100u64, 200, 300, 999] {
            let km = CanonicalKmer::new(raw);
            assert_eq!(table2.key_exists(km), table.key_exists(km));
            for uid in 0..3 {
                assert_eq!(
                    table2.key_occurs_in_unitig(km, uid),
                    table.key_occurs_in_unitig(km, uid)
                );
            }
        }

        // Verify specific positional queries
        let km300 = CanonicalKmer::new(300);
        assert!(table2.key_occurs_in_unitig_between(km300, 0, 20, 30));
        assert!(!table2.key_occurs_in_unitig_between(km300, 0, 21, 24));
    }

    #[test]
    fn test_serialization_empty() {
        let table = PoisonTable::build_from_occs(Vec::new()).unwrap();

        let mut buf = Vec::new();
        table.save(&mut buf).unwrap();

        let table2 = PoisonTable::load(&mut &buf[..]).unwrap();
        assert!(table2.empty());
        assert_eq!(table2.num_poison_kmers(), 0);
    }

    #[test]
    fn test_invalid_magic() {
        let data = b"BADMAGIC\x00\x00\x00\x00\x00\x00\x00\x00";
        let result = PoisonTable::load(&mut &data[..]);
        assert!(result.is_err());
        let err = result.unwrap_err().to_string();
        assert!(err.contains("invalid poison table magic"));
    }

    #[test]
    fn test_stats_json() {
        let table = PoisonTable::build_from_occs(make_test_occs()).unwrap();

        let mut buf = Vec::new();
        table.save_stats_json(&mut buf).unwrap();

        let json: serde_json::Value = serde_json::from_slice(&buf).unwrap();
        assert_eq!(json["num_poison_kmers"], 3);
        assert_eq!(json["num_poison_occs"], 6);
        assert_eq!(json["max_poison_occ"], 3);
    }

    #[test]
    fn test_single_kmer_single_occ() {
        let occs = vec![LabeledPoisonOcc {
            canonical_kmer: CanonicalKmer::new(42),
            unitig_id: 7,
            unitig_pos: 15,
        }];

        let table = PoisonTable::build_from_occs(occs).unwrap();

        assert_eq!(table.num_poison_kmers(), 1);
        assert_eq!(table.num_poison_occs(), 1);
        let km42 = CanonicalKmer::new(42);
        assert!(table.key_exists(km42));
        assert!(table.key_occurs_in_unitig(km42, 7));
        assert!(table.key_occurs_in_unitig_between(km42, 7, 15, 15));
        assert!(!table.key_occurs_in_unitig_between(km42, 7, 0, 14));
    }

    #[test]
    fn test_build_unsorted_input() {
        // Input in reverse order — should still work after sorting
        let occs = vec![
            LabeledPoisonOcc {
                canonical_kmer: CanonicalKmer::new(300),
                unitig_id: 0,
                unitig_pos: 30,
            },
            LabeledPoisonOcc {
                canonical_kmer: CanonicalKmer::new(100),
                unitig_id: 1,
                unitig_pos: 10,
            },
            LabeledPoisonOcc {
                canonical_kmer: CanonicalKmer::new(200),
                unitig_id: 2,
                unitig_pos: 0,
            },
            LabeledPoisonOcc {
                canonical_kmer: CanonicalKmer::new(100),
                unitig_id: 0,
                unitig_pos: 5,
            },
        ];

        let table = PoisonTable::build_from_occs(occs).unwrap();

        assert_eq!(table.num_poison_kmers(), 3);
        assert_eq!(table.num_poison_occs(), 4);
        assert!(table.key_occurs_in_unitig(CanonicalKmer::new(100), 0));
        assert!(table.key_occurs_in_unitig(CanonicalKmer::new(100), 1));
        assert!(table.key_occurs_in_unitig(CanonicalKmer::new(200), 2));
        assert!(table.key_occurs_in_unitig(CanonicalKmer::new(300), 0));
    }
}
