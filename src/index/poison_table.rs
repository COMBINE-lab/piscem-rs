//! Poison k-mer table — improves mapping specificity by identifying k-mers
//! from decoy sequences that sit at the boundary between present/absent on
//! reference unitigs.
//!
//! During mapping, if a read's hits contain a poison k-mer, the mapping is
//! discarded. This avoids false-positive pseudoalignments without requiring
//! alignment scoring.
//!
//! The table maps canonical k-mer values to their occurrence positions on
//! unitigs, in a single open-addressing array (see [`PoisonTable`]).
//!
//! Corresponds to the C++ `poison_table` class in `poison_table.hpp`.

use anyhow::{Context, Result, bail, ensure};
use std::io::{BufWriter, Read, Write};
use std::sync::OnceLock;
use tracing::{info, warn};

use crate::mapping::kmer_value::CanonicalKmer;

// ---------------------------------------------------------------------------
// Read-side pre-filter
// ---------------------------------------------------------------------------

/// Reverse the order of the low `k` 2-bit groups of `x` (`k <= 32`).
#[inline]
pub(crate) fn reverse_groups(x: u64, k: u32) -> u64 {
    let r = x.reverse_bits();
    // Full bit reversal also swaps the two bits inside each group; undo that.
    let r = ((r & 0x5555_5555_5555_5555) << 1) | ((r >> 1) & 0x5555_5555_5555_5555);
    r >> (64 - 2 * k)
}

/// Sub-linear candidate enumeration for small poison tables.
///
/// The read filter costs a roll and a probe at every one of the read's
/// k-mer positions. For a small table — e.g. the 38 poison k-mers of 10x
/// Flex `included=FALSE` probes — that per-position work is most of the
/// remaining poison cost, and almost all of it answers "no".
///
/// Instead, sample the read's q-grams (read encoding, `q < k`) at every
/// `stride = k - q + 1` positions. Any k-mer window `[p, p + k)` contains a
/// sampled q-gram, at the sample `s` with `p <= s <= p + k - q` (an interval
/// of `stride` positions), and if the window is a poison k-mer (in either
/// orientation) that q-gram is the pattern's q-gram at offset `s - p`. So a
/// bitmap over all patterns' q-grams, plus the list of `(pattern, offset)`
/// pairs behind each one, enumerates every window that could be poison; each
/// is then confirmed by a word compare and the exact query. A 150 bp read
/// needs ~10 bitmap tests instead of 128 probes.
///
/// Only built when the bitmap is sparse — the worst-case q-gram density is at
/// most [`QgramSampler::MAX_DENSITY`] for some `q <= 12` (a bitmap of at most
/// 2 MiB) — i.e. for up to a few thousand poison k-mers.
pub(crate) struct QgramSampler {
    pub(crate) q: u32,
    pub(crate) stride: usize,
    /// Bit `c` set iff some pattern has q-gram code `c`.
    bitmap: Vec<u64>,
    /// `(q-gram code, offset in pattern, pattern word)`, sorted by code.
    entries: Vec<(u32, u32, u64)>,
}

impl QgramSampler {
    const MIN_Q: u32 = 8;
    const MAX_Q: u32 = 12;
    /// Worst-case fraction of bitmap bits set.
    const MAX_DENSITY: usize = 128;

    /// The q-gram length to sample with for `n_patterns` read-encoded
    /// patterns (both orientations of each poison k-mer), if any is sparse
    /// enough to pay off.
    fn q_for(n_patterns: usize, k: u32) -> Option<u32> {
        (Self::MIN_Q.min(k)..=Self::MAX_Q.min(k)).find(|&q| {
            let max_codes = n_patterns * (k - q + 1) as usize;
            max_codes * Self::MAX_DENSITY <= 1usize << (2 * q)
        })
    }

    /// Build a sampler with q-gram length `q` (`1 <= q <= k`).
    pub(crate) fn build(words: &[u64], k: u32, q: u32) -> Self {
        assert!(1 <= q && q <= k && q <= 16, "q-gram length {q} for k={k}");
        let qmask = (1u64 << (2 * q)) - 1;
        let mut entries: Vec<(u32, u32, u64)> = words
            .iter()
            .flat_map(|&w| (0..=k - q).map(move |off| (((w >> (2 * off)) & qmask) as u32, off, w)))
            .collect();
        entries.sort_unstable();
        entries.dedup();
        let mut bitmap = vec![0u64; ((1usize << (2 * q)) / 64).max(1)];
        for &(code, _, _) in &entries {
            bitmap[code as usize / 64] |= 1u64 << (code % 64);
        }
        Self {
            q,
            stride: (k - q + 1) as usize,
            bitmap,
            entries,
        }
    }

    /// Whether some pattern has q-gram `code`.
    #[inline(always)]
    pub(crate) fn has(&self, code: u32) -> bool {
        // SAFETY: `code < 4^q` and the bitmap holds `4^q` bits.
        let w = unsafe { *self.bitmap.get_unchecked(code as usize / 64) };
        (w >> (code % 64)) & 1 != 0
    }

    /// The `(offset, pattern word)` pairs behind q-gram `code`.
    #[inline]
    pub(crate) fn patterns_with(&self, code: u32) -> impl Iterator<Item = (u32, u64)> + '_ {
        let start = self.entries.partition_point(|e| e.0 < code);
        self.entries[start..]
            .iter()
            .take_while(move |e| e.0 == code)
            .map(|e| (e.1, e.2))
    }
}

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
    /// Sampled-scan alternative for small tables; see [`QgramSampler`].
    pub(crate) sampler: Option<QgramSampler>,
    /// `64 - log2(words.len())`: shift that maps a hash to a word index.
    word_shift: u32,
    k: u32,
}

impl ReadKmerFilter {
    const MIN_LOG2_WORDS: u32 = 9;
    const BITS_PER_ENTRY: usize = 16;

    /// Build the filter (and, for small tables, the sampler) over the
    /// `n_keys` canonical poison k-mers produced by `keys`.
    fn build(keys: impl Iterator<Item = CanonicalKmer> + Clone, n_keys: usize, k: u32) -> Self {
        assert!(
            (1..=31).contains(&k),
            "poison k-mers support 1 <= k <= 31, got {k}"
        );
        let words = |keys: &mut dyn Iterator<Item = CanonicalKmer>| -> Vec<u64> {
            keys.flat_map(|km| {
                let (fw, rc) = Self::read_words_of(km, k);
                [fw, rc]
            })
            .collect()
        };
        let n_patterns = 2 * n_keys;
        // Only small tables get a sampler, so only they materialize the
        // pattern list; a genome-scale table streams straight into the filter.
        let sampler = QgramSampler::q_for(n_patterns, k)
            .map(|q| QgramSampler::build(&words(&mut keys.clone()), k, q));
        let n_words = (n_patterns * Self::BITS_PER_ENTRY / 64)
            .next_power_of_two()
            .max(1 << Self::MIN_LOG2_WORDS);
        let mut f = Self {
            words: vec![0u64; n_words],
            sampler,
            word_shift: 64 - n_words.trailing_zeros(),
            k,
        };
        for km in keys {
            let (fw, rc) = Self::read_words_of(km, k);
            for word in [fw, rc] {
                let (w, bits) = f.slot(word);
                f.words[w] |= bits;
            }
        }
        match &f.sampler {
            Some(s) => info!(
                "poison scan: sampled ({n_patterns} patterns, q={}, stride {})",
                s.q, s.stride
            ),
            None => info!(
                "poison scan: every position ({n_patterns} patterns, {} KiB filter)",
                n_words * 8 / 1024
            ),
        }
        f
    }

    /// The two read-encoding words a canonical k-mer can appear as.
    ///
    /// `km` is packed first-base-most-significant with `A=0 C=1 G=2 T=3`.
    /// Register-only, so building the filter for millions of keys is cheap.
    pub(crate) fn read_words_of(km: CanonicalKmer, k: u32) -> (u64, u64) {
        // ACGT code -> read code (A=0 C=1 T=2 G=3) in every 2-bit group: swap
        // G and T by flipping the low bit wherever the high bit is set.
        let v = km.as_u64();
        let v = v ^ ((v >> 1) & 0x5555_5555_5555_5555);
        // Read as LSB-first, `v` spells the k-mer backwards; reversing the
        // group order gives the forward word. Backwards *and* complemented
        // (read-encoding complement is `^ 0b10`) is the reverse complement.
        let fw = reverse_groups(v, k);
        let rc = v ^ (0xAAAA_AAAA_AAAA_AAAA >> (64 - 2 * k));
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

impl PoisonOcc {
    fn from_le_bytes(b: &[u8; 8]) -> Self {
        Self {
            unitig_id: u32::from_le_bytes(b[..4].try_into().unwrap()),
            unitig_pos: u32::from_le_bytes(b[4..].try_into().unwrap()),
        }
    }
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

/// Magic bytes of the original piscem-rs poison table format (offsets,
/// occurrences, then sorted `(k-mer, offset index)` map entries). Still read.
const POISON_MAGIC_V1: &[u8; 8] = b"PPOIS01\0";

/// Magic bytes of the flat poison table format: the [`PoisonTable`] slot
/// array and overflow occurrences verbatim, so loading is a bulk read.
const POISON_MAGIC_V2: &[u8; 8] = b"PPOIS02\0";

/// Seed of layout attempt `attempt`; see [`PoisonTable::layout`]. Every
/// writer tries them in this order, so a key set always gets the same seed
/// (almost always attempt 0's) and the same file.
const fn seed_for(attempt: u64) -> u64 {
    (attempt + 1).wrapping_mul(0x9E37_79B9_7F4A_7C15)
}

/// Seeds tried before placing without a probe budget. Never reached in
/// practice: no seed of 6,200 tried on real decoy k-mers (40 to 3.7M keys)
/// exceeded the budget, and a seed only fails on a key set built against it.
const MAX_LAYOUT_ATTEMPTS: u64 = 16;

/// murmur3's 64-bit finalizer: a bijection with full avalanche.
#[inline(always)]
fn fmix64(mut k: u64) -> u64 {
    k ^= k >> 33;
    k = k.wrapping_mul(0xff51_afd7_ed55_8ccd);
    k ^= k >> 33;
    k = k.wrapping_mul(0xc4ce_b9fe_1a85_ec53);
    k ^ (k >> 33)
}

/// Home slot of `key` in a table of `num_slots` slots under `seed`:
/// fastrange over `fmix64(key ^ seed)`.
#[inline(always)]
fn home_slot(key: u64, seed: u64, num_slots: usize) -> usize {
    ((u128::from(fmix64(key ^ seed)) * num_slots as u128) >> 64) as usize
}

/// How much probing a layout may need before it is rejected.
///
/// With a random hash, linear probing at load 2/3 displaces a key by 1 slot on
/// average (2 probes per hit, 5 per miss) and the longest displacement for
/// 3.7M keys is ~100-150. The budget is 3x the mean plus slack for small
/// tables, and a longest displacement of `64 + 16 log2(n)` (~410 at 3.7M
/// keys), so a layout within it queries like a random hash's and a seed
/// that clusters some keys is caught.
struct ProbeBudget {
    total: u64,
    max: usize,
}

impl ProbeBudget {
    fn for_keys(n: usize) -> Self {
        Self {
            total: 3 * n as u64 + 1024,
            max: 64 + 16 * (usize::BITS - n.leading_zeros()) as usize,
        }
    }
}

/// One slot of the [`PoisonTable`] array.
///
/// `key` is the canonical k-mer (< 2^62), with [`Slot::MULTI`] set when the
/// k-mer has more than one occurrence, or [`Slot::EMPTY`]. A single
/// occurrence — 96% of keys for a genome decoy — is stored inline
/// (`a` = unitig, `b` = position); otherwise `a..a + b` indexes the overflow
/// occurrences.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Slot {
    key: u64,
    a: u32,
    b: u32,
}

impl Slot {
    const EMPTY: Slot = Slot {
        key: u64::MAX,
        a: 0,
        b: 0,
    };
    const MULTI: u64 = 1 << 63;
    const BYTES: usize = 16;

    #[inline(always)]
    fn kmer(&self) -> u64 {
        self.key & !Self::MULTI
    }

    fn to_le_bytes(self) -> [u8; Self::BYTES] {
        let mut out = [0u8; Self::BYTES];
        out[..8].copy_from_slice(&self.key.to_le_bytes());
        out[8..12].copy_from_slice(&self.a.to_le_bytes());
        out[12..].copy_from_slice(&self.b.to_le_bytes());
        out
    }

    fn from_le_bytes(b: &[u8; Self::BYTES]) -> Self {
        Self {
            key: u64::from_le_bytes(b[..8].try_into().unwrap()),
            a: u32::from_le_bytes(b[8..12].try_into().unwrap()),
            b: u32::from_le_bytes(b[12..16].try_into().unwrap()),
        }
    }
}

/// Poison k-mer table for improving mapping specificity.
///
/// Maps canonical k-mer values to their occurrence positions on reference
/// unitigs. Lookups only happen for read k-mers that pass the read-side
/// pre-filter, so each should cost as few cache misses as possible: keys live
/// in one open-addressing array (linear probing, load factor <= 2/3,
/// fastrange slot selection), with a single occurrence stored in the slot
/// itself — so a typical query is one cache line, against four dependent
/// accesses (hash control bytes, map slot, offsets, occurrences) for a hash
/// map of offset indices. Slot placement depends only on the keys (inserted
/// in sorted order) and the hash seed, so the array is deterministic and
/// serialized as is.
///
/// The seed is chosen per key set, not fixed: [`PoisonTable::layout`] keeps
/// the first seed whose layout probes no more than a random hash would (see
/// [`ProbeBudget`]), and loading re-checks that bound, rebuilding with a new
/// seed if a file does not meet it. So no key set, however structured, can
/// leave the table slow to build, load or query.
pub struct PoisonTable {
    slots: Vec<Slot>,
    /// Home-slot hash seed; see [`home_slot`].
    seed: u64,
    /// Occurrences of the k-mers with more than one, grouped by k-mer.
    overflow: Vec<PoisonOcc>,
    num_kmers: usize,
    num_occs: usize,
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
    /// 3. Group the occurrences by k-mer
    pub fn build_from_occs(mut occs: Vec<LabeledPoisonOcc>) -> Result<Self> {
        occs.sort_unstable_by(|a, b| {
            a.canonical_kmer
                .cmp(&b.canonical_kmer)
                .then(a.unitig_id.cmp(&b.unitig_id))
                .then(a.unitig_pos.cmp(&b.unitig_pos))
        });
        occs.dedup();

        info!("Total number of distinct poison k-mer occs: {}", occs.len());

        let groups: Vec<&[LabeledPoisonOcc]> = occs
            .chunk_by(|a, b| a.canonical_kmer == b.canonical_kmer)
            .collect();
        let table = Self::from_groups(
            groups.iter().map(|g| {
                (
                    g[0].canonical_kmer,
                    g.iter().map(|o| PoisonOcc {
                        unitig_id: o.unitig_id,
                        unitig_pos: o.unitig_pos,
                    }),
                )
            }),
            occs.len(),
        )?;
        info!(
            "[poison_table]: The most frequently occurring poison k-mer appeared in {} distinct unitig positions.",
            table.max_poison_occ()
        );
        Ok(table)
    }

    /// Build from occurrence groups, one per distinct k-mer, in increasing
    /// k-mer order (which fixes the slot layout).
    fn from_groups<I, O>(groups: I, num_occs: usize) -> Result<Self>
    where
        I: ExactSizeIterator<Item = (CanonicalKmer, O)>,
        O: Iterator<Item = PoisonOcc>,
    {
        ensure!(
            num_occs <= u32::MAX as usize,
            "Too many poison k-mer occurrences ({num_occs}) for u32 overflow indices"
        );
        let num_kmers = groups.len();
        let mut entries = Vec::with_capacity(num_kmers);
        let mut overflow = Vec::new();
        let mut prev: Option<u64> = None;
        for (km, mut occs) in groups {
            let key = km.as_u64();
            ensure!(key < 1 << 62, "poison k-mer value {key:#x} exceeds 62 bits");
            ensure!(
                prev.is_none_or(|p| p < key),
                "poison k-mer groups must be distinct and sorted"
            );
            prev = Some(key);
            let first = occs.next().context("poison k-mer without occurrences")?;
            let slot = match occs.next() {
                None => Slot {
                    key,
                    a: first.unitig_id,
                    b: first.unitig_pos,
                },
                Some(second) => {
                    let start = overflow.len();
                    overflow.extend([first, second]);
                    overflow.extend(occs);
                    Slot {
                        key: key | Slot::MULTI,
                        a: start as u32,
                        b: (overflow.len() - start) as u32,
                    }
                }
            };
            entries.push(slot);
        }
        let inline = entries.iter().filter(|s| s.key & Slot::MULTI == 0).count();
        ensure!(
            inline + overflow.len() == num_occs,
            "poison table holds {} occurrences, expected {num_occs}",
            inline + overflow.len()
        );
        let (slots, seed) = Self::layout(&entries);
        Ok(Self {
            slots,
            seed,
            overflow,
            num_kmers,
            num_occs,
            read_filter: OnceLock::new(),
        })
    }

    /// Place `entries` (distinct keys, in increasing k-mer order) under the
    /// first seed of [`seed_for`]'s sequence whose layout stays within the
    /// [`ProbeBudget`]. An attempt is abandoned as soon as it exceeds the
    /// budget, so even a seed that sends every key to one slot costs O(n).
    fn layout(entries: &[Slot]) -> (Vec<Slot>, u64) {
        let budget = ProbeBudget::for_keys(entries.len());
        for attempt in 0..MAX_LAYOUT_ATTEMPTS {
            let seed = seed_for(attempt);
            if let Some(slots) = Self::place(entries, seed, Some(&budget)) {
                if attempt > 0 {
                    info!(
                        "[poison_table]: hash seed {attempt} of the sequence used; \
                         earlier seeds clustered these {} k-mers",
                        entries.len()
                    );
                }
                return (slots, seed);
            }
        }
        let seed = seed_for(MAX_LAYOUT_ATTEMPTS);
        let slots = Self::place(entries, seed, None).expect("unbudgeted placement cannot fail");
        (slots, seed)
    }

    /// Linear-probing insertion of `entries` in order, or `None` once the
    /// layout exceeds `budget`.
    fn place(entries: &[Slot], seed: u64, budget: Option<&ProbeBudget>) -> Option<Vec<Slot>> {
        let cap = Self::capacity_for(entries.len());
        let mut slots = vec![Slot::EMPTY; cap];
        let mut total = 0u64;
        for e in entries {
            let mut i = home_slot(e.kmer(), seed, cap);
            let mut d = 0usize;
            while slots[i] != Slot::EMPTY {
                i = if i + 1 == cap { 0 } else { i + 1 };
                d += 1;
                if budget.is_some_and(|b| d > b.max) {
                    return None;
                }
            }
            total += d as u64;
            if budget.is_some_and(|b| total > b.total) {
                return None;
            }
            slots[i] = *e;
        }
        Some(slots)
    }

    /// Check a loaded slot array: every key reachable from its home slot (else
    /// the file is corrupt) and every overflow range in bounds. Returns
    /// whether the layout is within the [`ProbeBudget`]. One linear pass, so
    /// it costs O(slots) however badly the keys are placed.
    fn check_layout(&self) -> Result<bool> {
        let cap = self.slots.len();
        let budget = ProbeBudget::for_keys(self.num_kmers);
        let start = self
            .slots
            .iter()
            .position(|s| *s == Slot::EMPTY)
            .context("corrupt poison table: no empty slot")?;
        // `run`: occupied slots immediately before the current one. A key is
        // reachable iff every slot from its home up to its own is occupied,
        // i.e. iff its displacement is at most `run`.
        let (mut run, mut total, mut max, mut n, mut inline) =
            (0usize, 0u64, 0usize, 0usize, 0usize);
        for j in 1..=cap {
            let p = (start + j) % cap;
            let s = &self.slots[p];
            if *s == Slot::EMPTY {
                run = 0;
                continue;
            }
            let km = s.kmer();
            ensure!(
                km < 1 << 62,
                "corrupt poison table: k-mer {km:#x} exceeds 62 bits"
            );
            let d = (p + cap - home_slot(km, self.seed, cap)) % cap;
            ensure!(
                d <= run,
                "corrupt poison table: k-mer {km:#x} not reachable from its home slot"
            );
            if s.key & Slot::MULTI == 0 {
                inline += 1;
            } else {
                ensure!(
                    s.b >= 2 && (s.a as usize + s.b as usize) <= self.overflow.len(),
                    "corrupt poison table: overflow range out of bounds"
                );
            }
            total += d as u64;
            max = max.max(d);
            n += 1;
            run += 1;
        }
        ensure!(
            n == self.num_kmers,
            "poison table holds {n} k-mers, header says {}",
            self.num_kmers
        );
        ensure!(
            inline + self.overflow.len() == self.num_occs,
            "poison table holds {} occurrences, header says {}",
            inline + self.overflow.len(),
            self.num_occs
        );
        Ok(total <= budget.total && max <= budget.max)
    }

    /// Slot count for `n` keys: load factor at most 2/3, never full.
    fn capacity_for(n: usize) -> usize {
        (n + n / 2).max(n + 1).max(8)
    }

    /// Home slot of `key` in this table.
    #[inline(always)]
    fn home(&self, key: u64) -> usize {
        home_slot(key, self.seed, self.slots.len())
    }

    #[inline]
    fn find(&self, km: CanonicalKmer) -> Option<&Slot> {
        if self.num_kmers == 0 {
            return None;
        }
        let key = km.as_u64();
        let cap = self.slots.len();
        let mut i = self.home(key);
        loop {
            let s = &self.slots[i];
            if s.key == Slot::EMPTY.key {
                return None;
            }
            if s.kmer() == key {
                return Some(s);
            }
            i = if i + 1 == cap { 0 } else { i + 1 };
        }
    }

    /// Whether any occurrence of the k-mer in `slot` satisfies `pred`.
    #[inline]
    fn any_occ(&self, slot: &Slot, pred: impl Fn(&PoisonOcc) -> bool) -> bool {
        if slot.key & Slot::MULTI == 0 {
            pred(&PoisonOcc {
                unitig_id: slot.a,
                unitig_pos: slot.b,
            })
        } else {
            self.overflow[slot.a as usize..(slot.a + slot.b) as usize]
                .iter()
                .any(pred)
        }
    }

    // -----------------------------------------------------------------------
    // Query methods
    // -----------------------------------------------------------------------

    /// Whether this poison table is empty (contains no poison k-mers).
    #[inline]
    pub fn empty(&self) -> bool {
        self.num_kmers == 0
    }

    /// Number of distinct poison k-mers.
    #[inline]
    pub fn num_poison_kmers(&self) -> usize {
        self.num_kmers
    }

    /// Total number of poison k-mer occurrences.
    #[inline]
    pub fn num_poison_occs(&self) -> usize {
        self.num_occs
    }

    /// Returns `true` if `km` is a poison k-mer.
    #[inline]
    pub fn key_exists(&self, km: CanonicalKmer) -> bool {
        self.find(km).is_some()
    }

    /// Returns `true` if poison k-mer `km` occurs on unitig `unitig_id`.
    #[inline]
    pub fn key_occurs_in_unitig(&self, km: CanonicalKmer, unitig_id: u32) -> bool {
        self.find(km)
            .is_some_and(|s| self.any_occ(s, |o| o.unitig_id == unitig_id))
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
        self.find(km).is_some_and(|s| {
            self.any_occ(s, |o| {
                o.unitig_id == unitig_id && o.unitig_pos >= lb && o.unitig_pos <= ub
            })
        })
    }

    /// The canonical poison k-mers, in slot order.
    pub fn keys(&self) -> impl Iterator<Item = CanonicalKmer> + Clone + '_ {
        self.slots
            .iter()
            .filter(|s| **s != Slot::EMPTY)
            .map(|s| CanonicalKmer::new(s.kmer()))
    }

    /// Maximum number of occurrences for any single poison k-mer.
    pub fn max_poison_occ(&self) -> u32 {
        self.slots
            .iter()
            .filter(|s| **s != Slot::EMPTY)
            .map(|s| if s.key & Slot::MULTI == 0 { 1 } else { s.b })
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
            .get_or_init(|| ReadKmerFilter::build(self.keys(), self.num_kmers, k));
        debug_assert_eq!(f.k, k, "poison read filter requested for a different k");
        f
    }

    // -----------------------------------------------------------------------
    // Serialization
    // -----------------------------------------------------------------------

    /// Serialize the poison table in the flat `PPOIS02` layout (the
    /// default since piscem-rs 0.11.0, whose indexes must be rebuilt anyway).
    /// It saves the load-time table rebuild and ~13% of the file size against
    /// `PPOIS01`, which [`PoisonTable::save_v1`] still writes. Releases before
    /// 0.11.0 cannot read it.
    ///
    /// Format (all little-endian):
    /// ```text
    /// [magic: 8 bytes "PPOIS02\0"]
    /// [num_kmers: u64] [num_occs: u64] [num_slots: u64] [num_overflow: u64]
    /// [seed: u64]
    /// [slots: num_slots × {key: u64, a: u32, b: u32}]
    /// [overflow: num_overflow × {unitig_id: u32, unitig_pos: u32}]
    /// ```
    /// Slot positions follow from [`home_slot`] under `seed`; the hash
    /// function is part of the format, the seed is per file.
    pub fn save<W: Write>(&self, writer: &mut W) -> Result<()> {
        let mut w = BufWriter::new(writer);
        w.write_all(POISON_MAGIC_V2)?;
        for n in [
            self.num_kmers,
            self.num_occs,
            self.slots.len(),
            self.overflow.len(),
        ] {
            w.write_all(&(n as u64).to_le_bytes())?;
        }
        w.write_all(&self.seed.to_le_bytes())?;
        for s in &self.slots {
            w.write_all(&s.to_le_bytes())?;
        }
        for o in &self.overflow {
            w.write_all(&o.unitig_id.to_le_bytes())?;
            w.write_all(&o.unitig_pos.to_le_bytes())?;
        }
        w.flush()?;
        Ok(())
    }

    /// Serialize the poison table in the original `PPOIS01` layout, which
    /// every piscem-rs release reads (and C++ piscem writes).
    /// [`PoisonTable::load`] converts it into the flat in-memory table;
    /// [`PoisonTable::save`] writes that table directly.
    ///
    /// ```text
    /// [magic: 8 bytes "PPOIS01\0"]
    /// [num_kmers: u64 LE]
    /// [num_offsets: u64 LE]
    /// [offsets: num_offsets × u32 LE]
    /// [num_occs: u64 LE]
    /// [poison_occs: num_occs × {unitig_id: u32 LE, unitig_pos: u32 LE}]
    /// [map_entries: num_kmers × {canonical_kmer: u64 LE, offset_idx: u64 LE}]
    /// ```
    pub fn save_v1<W: Write>(&self, writer: &mut W) -> Result<()> {
        let mut keyed: Vec<&Slot> = self.slots.iter().filter(|s| **s != Slot::EMPTY).collect();
        keyed.sort_unstable_by_key(|s| s.kmer());
        let mut w = BufWriter::new(writer);
        w.write_all(POISON_MAGIC_V1)?;
        w.write_all(&(self.num_kmers as u64).to_le_bytes())?;
        w.write_all(&(self.num_kmers as u64 + 1).to_le_bytes())?;
        let mut off = 0u32;
        w.write_all(&off.to_le_bytes())?;
        for s in &keyed {
            off += if s.key & Slot::MULTI == 0 { 1 } else { s.b };
            w.write_all(&off.to_le_bytes())?;
        }
        w.write_all(&(self.num_occs as u64).to_le_bytes())?;
        for s in &keyed {
            let mut write = |o: &PoisonOcc| -> std::io::Result<()> {
                w.write_all(&o.unitig_id.to_le_bytes())?;
                w.write_all(&o.unitig_pos.to_le_bytes())
            };
            if s.key & Slot::MULTI == 0 {
                write(&PoisonOcc {
                    unitig_id: s.a,
                    unitig_pos: s.b,
                })?;
            } else {
                for o in &self.overflow[s.a as usize..(s.a + s.b) as usize] {
                    write(o)?;
                }
            }
        }
        for (idx, s) in keyed.iter().enumerate() {
            w.write_all(&s.kmer().to_le_bytes())?;
            w.write_all(&(idx as u64).to_le_bytes())?;
        }
        w.flush()?;
        Ok(())
    }

    /// Deserialize a poison table in either format (`PPOIS02`, or the
    /// original `PPOIS01`, which is converted).
    pub fn load<R: Read>(reader: &mut R) -> Result<Self> {
        let mut magic = [0u8; 8];
        reader
            .read_exact(&mut magic)
            .context("failed to read poison table magic")?;
        if magic == *POISON_MAGIC_V2 {
            Self::load_v2(reader)
        } else if magic == *POISON_MAGIC_V1 {
            Self::load_v1(reader)
        } else {
            bail!(
                "invalid poison table magic: expected {:?} or {:?}, got {:?}",
                POISON_MAGIC_V2,
                POISON_MAGIC_V1,
                magic
            );
        }
    }

    fn load_v2<R: Read>(r: &mut R) -> Result<Self> {
        let mut header = [0u64; 5];
        for h in &mut header {
            *h = read_u64_le(r).context("failed to read poison table header")?;
        }
        let seed = header[4];
        let [num_kmers, num_occs, num_slots, num_overflow] =
            [header[0], header[1], header[2], header[3]].map(|h| h as usize);
        ensure!(
            num_slots == Self::capacity_for(num_kmers) && num_overflow <= num_occs,
            "inconsistent poison table header ({num_kmers} k-mers, {num_occs} occurrences, \
             {num_slots} slots, {num_overflow} overflow)"
        );
        let slots = read_records(r, num_slots, Slot::from_le_bytes)
            .context("failed to read poison slots")?;
        let overflow = read_records(r, num_overflow, PoisonOcc::from_le_bytes)
            .context("failed to read poison overflow occurrences")?;
        let mut table = Self {
            slots,
            seed,
            overflow,
            num_kmers,
            num_occs,
            read_filter: OnceLock::new(),
        };
        // Every key must be findable from its home slot, and every overflow
        // range in bounds, or queries could miss or panic.
        if table.check_layout()? {
            // Within the probe budget, so this pass is linear too. It catches
            // duplicate keys, which `find` would resolve to the first copy.
            for s in table.slots.iter().filter(|s| **s != Slot::EMPTY) {
                ensure!(
                    table.find(CanonicalKmer::new(s.kmer())) == Some(s),
                    "corrupt poison table: duplicate k-mer {:#x}",
                    s.kmer()
                );
            }
        } else {
            // Valid but clustered (not something `save` writes): lay the
            // keys out again rather than query a slow table.
            let mut entries: Vec<Slot> = table
                .slots
                .iter()
                .copied()
                .filter(|s| *s != Slot::EMPTY)
                .collect();
            entries.sort_unstable_by_key(|s| s.kmer());
            ensure!(
                entries.windows(2).all(|w| w[0].kmer() < w[1].kmer()),
                "corrupt poison table: duplicate k-mers"
            );
            let (slots, seed) = Self::layout(&entries);
            warn!(
                "poison table layout (seed {:#x}) exceeds the probe budget; re-laid out in memory",
                table.seed
            );
            table.slots = slots;
            table.seed = seed;
        }
        Ok(table)
    }

    fn load_v1<R: Read>(r: &mut R) -> Result<Self> {
        let num_kmers = read_u64_le(r).context("failed to read num_kmers")? as usize;
        let num_offsets = read_u64_le(r).context("failed to read num_offsets")? as usize;
        let offsets = read_records(r, num_offsets, |b: &[u8; 4]| u32::from_le_bytes(*b))
            .context("failed to read poison offsets")?;
        let num_occs = read_u64_le(r).context("failed to read num_occs")? as usize;
        let occs = read_records(r, num_occs, PoisonOcc::from_le_bytes)
            .context("failed to read poison occurrences")?;
        // `(k-mer, offset index)` entries -> `(k-mer, occurrence range)`.
        let mut groups = read_records(r, num_kmers, |e: &[u8; 16]| {
            let km = u64::from_le_bytes(e[..8].try_into().unwrap());
            let idx = u64::from_le_bytes(e[8..].try_into().unwrap());
            (km, idx)
        })
        .context("failed to read poison map entries")?
        .into_iter()
        .map(|(km, idx)| {
            let idx = usize::try_from(idx).ok().filter(|&i| i + 1 < offsets.len());
            let idx = idx.context("poison map index out of range")?;
            let (lo, hi) = (offsets[idx], offsets[idx + 1]);
            ensure!(
                lo < hi && hi as usize <= num_occs,
                "bad poison occurrence range"
            );
            Ok((km, lo, hi))
        })
        .collect::<Result<Vec<_>>>()?;
        drop(offsets);
        // Every writer emits entries in k-mer order; sort only if not.
        if !groups.is_sorted_by_key(|g| g.0) {
            groups.sort_unstable_by_key(|g| g.0);
        }
        Self::from_groups(
            groups.iter().map(|&(km, lo, hi)| {
                (
                    CanonicalKmer::new(km),
                    occs[lo as usize..hi as usize].iter().copied(),
                )
            }),
            groups.iter().map(|g| (g.2 - g.1) as usize).sum(),
        )
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

/// Read `count` fixed-size little-endian records through a bounded buffer,
/// so the raw bytes never coexist in full with the decoded records.
fn read_records<R: Read, T, const N: usize>(
    reader: &mut R,
    count: usize,
    decode: impl Fn(&[u8; N]) -> T,
) -> std::io::Result<Vec<T>> {
    const CHUNK_BYTES: usize = 1 << 20;
    let per_chunk = (CHUNK_BYTES / N).max(1);
    let mut out = Vec::with_capacity(count);
    let mut buf = vec![0u8; per_chunk.min(count) * N];
    let mut left = count;
    while left > 0 {
        let n = left.min(per_chunk);
        reader.read_exact(&mut buf[..n * N])?;
        out.extend(buf[..n * N].as_chunks::<N>().0.iter().map(&decode));
        left -= n;
    }
    Ok(out)
}

fn read_u64_le<R: Read>(reader: &mut R) -> std::io::Result<u64> {
    let mut buf = [0u8; 8];
    reader.read_exact(&mut buf)?;
    Ok(u64::from_le_bytes(buf))
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Brute-force model of the table queries over a flat occurrence list.
    struct Model(Vec<LabeledPoisonOcc>);
    impl Model {
        fn exists(&self, km: u64) -> bool {
            self.0.iter().any(|o| o.canonical_kmer.as_u64() == km)
        }
        fn on(&self, km: u64, u: u32) -> bool {
            self.0
                .iter()
                .any(|o| o.canonical_kmer.as_u64() == km && o.unitig_id == u)
        }
        fn between(&self, km: u64, u: u32, lb: u32, ub: u32) -> bool {
            self.0.iter().any(|o| {
                o.canonical_kmer.as_u64() == km
                    && o.unitig_id == u
                    && o.unitig_pos >= lb
                    && o.unitig_pos <= ub
            })
        }
    }

    fn assert_matches_model(t: &PoisonTable, m: &Model, keyspace: u64, what: &str) {
        for km in 0..keyspace {
            let c = CanonicalKmer::new(km);
            assert_eq!(t.key_exists(c), m.exists(km), "{what}: exists {km}");
            for u in 0..4 {
                assert_eq!(
                    t.key_occurs_in_unitig(c, u),
                    m.on(km, u),
                    "{what}: on {km} {u}"
                );
                for (lb, ub) in [(0, u32::MAX), (3, 9), (5, 5), (10, 2), (0, 0)] {
                    assert_eq!(
                        t.key_occurs_in_unitig_between(c, u, lb, ub),
                        m.between(km, u, lb, ub),
                        "{what}: between {km} {u} {lb} {ub}"
                    );
                }
            }
        }
    }

    #[test]
    fn test_flat_table_matches_model_and_round_trips() {
        let mut x = 0x0123_4567_89ab_cdefu64;
        let mut next = move || {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            x
        };
        for n in [0usize, 1, 2, 7, 50, 400] {
            // Keys from a small space (so many repeat, giving multi-occurrence
            // k-mers) plus some large values near the 62-bit limit.
            let occs: Vec<LabeledPoisonOcc> = (0..n)
                .map(|i| LabeledPoisonOcc {
                    canonical_kmer: CanonicalKmer::new(if i % 9 == 0 {
                        (1u64 << 62) - 1 - (next() % 3)
                    } else {
                        next() % 200
                    }),
                    unitig_id: (next() % 4) as u32,
                    unitig_pos: (next() % 12) as u32,
                })
                .collect();
            let model = Model(occs.clone());
            let t = PoisonTable::build_from_occs(occs).unwrap();
            let mut distinct: Vec<u64> =
                model.0.iter().map(|o| o.canonical_kmer.as_u64()).collect();
            distinct.sort_unstable();
            distinct.dedup();
            assert_eq!(t.num_poison_kmers(), distinct.len());
            let mut keys: Vec<u64> = t.keys().map(|k| k.as_u64()).collect();
            keys.sort_unstable();
            assert_eq!(keys, distinct);
            assert_matches_model(&t, &model, 210, &format!("built n={n}"));

            let mut v2 = Vec::new();
            t.save(&mut v2).unwrap();
            let t2 = PoisonTable::load(&mut v2.as_slice()).unwrap();
            assert_matches_model(&t2, &model, 210, &format!("PPOIS02 n={n}"));
            assert_eq!(t2.max_poison_occ(), t.max_poison_occ());
            assert_eq!(t2.num_poison_occs(), t.num_poison_occs());

            let mut v1 = Vec::new();
            t.save_v1(&mut v1).unwrap();
            let t1 = PoisonTable::load(&mut v1.as_slice()).unwrap();
            assert_matches_model(&t1, &model, 210, &format!("PPOIS01 n={n}"));
            // Saving is deterministic in both formats.
            let (mut v2b, mut v1b) = (Vec::new(), Vec::new());
            t1.save(&mut v2b).unwrap();
            t2.save_v1(&mut v1b).unwrap();
            assert_eq!(v2b, v2, "PPOIS02 bytes n={n}");
            assert_eq!(v1b, v1, "PPOIS01 bytes n={n}");
        }
    }

    #[test]
    fn test_default_save_writes_ppois02() {
        let occs = vec![LabeledPoisonOcc {
            canonical_kmer: CanonicalKmer::new(42),
            unitig_id: 3,
            unitig_pos: 9,
        }];
        let t = PoisonTable::build_from_occs(occs).unwrap();
        let (mut v2, mut v1) = (Vec::new(), Vec::new());
        t.save(&mut v2).unwrap();
        t.save_v1(&mut v1).unwrap();
        assert_eq!(&v2[..8], POISON_MAGIC_V2);
        assert_eq!(&v1[..8], POISON_MAGIC_V1);
    }

    /// Inverse of [`fmix64`], to build key sets a given seed sends to one slot.
    fn unfmix64(mut k: u64) -> u64 {
        let inv = |a: u64| {
            let mut x = a;
            for _ in 0..6 {
                x = x.wrapping_mul(2u64.wrapping_sub(a.wrapping_mul(x)));
            }
            x
        };
        k ^= k >> 33;
        k = k.wrapping_mul(inv(0xc4ce_b9fe_1a85_ec53));
        k ^= k >> 33;
        k = k.wrapping_mul(inv(0xff51_afd7_ed55_8ccd));
        k ^ (k >> 33)
    }

    /// `n` distinct 62-bit keys that all hash to slot 0 under `seed`.
    fn colliding_keys(seed: u64, n: usize) -> Vec<u64> {
        let keys: Vec<u64> = (1u64..)
            .map(|i| unfmix64(i) ^ seed)
            .filter(|&k| k < 1 << 62)
            .take(n)
            .collect();
        for &k in &keys {
            assert_eq!(fmix64(k ^ seed) >> 40, 0);
        }
        keys
    }

    /// Queries on a table built from `occs_of(keys)`, plus absent keys that
    /// collide with them.
    fn check_occs(t: &PoisonTable, keys: &[u64], absent: &[u64], what: &str) {
        for (i, &k) in keys.iter().enumerate() {
            let c = CanonicalKmer::new(k);
            let u = i as u32;
            assert!(t.key_exists(c), "{what}: key {i}");
            assert!(t.key_occurs_in_unitig(c, u), "{what}: unitig {i}");
            assert!(
                !t.key_occurs_in_unitig(c, u + 1),
                "{what}: other unitig {i}"
            );
            assert!(
                t.key_occurs_in_unitig_between(c, u, 0, 0),
                "{what}: pos 0 {i}"
            );
            assert_eq!(
                t.key_occurs_in_unitig_between(c, u, 1, 1),
                i % 3 == 0,
                "{what}: pos 1 {i}"
            );
        }
        for &k in absent {
            assert!(
                !t.key_exists(CanonicalKmer::new(k)),
                "{what}: absent {k:#x}"
            );
        }
    }

    fn occs_of(keys: &[u64]) -> Vec<LabeledPoisonOcc> {
        keys.iter()
            .enumerate()
            .flat_map(|(i, &k)| {
                // every third k-mer gets a second occurrence
                (0..1 + (i % 3 == 0) as u32).map(move |j| LabeledPoisonOcc {
                    canonical_kmer: CanonicalKmer::new(k),
                    unitig_id: i as u32,
                    unitig_pos: j,
                })
            })
            .collect()
    }

    /// A key set the first seed piles onto one slot (quadratic to build and
    /// to query) is laid out under the next seed, deterministically, and
    /// within the probe budget.
    #[test]
    fn test_clustering_key_set_gets_another_seed() {
        let all = colliding_keys(seed_for(0), 5100);
        let (keys, absent) = all.split_at(5000);
        let occs = occs_of(keys);
        let t = PoisonTable::build_from_occs(occs.clone()).unwrap();
        assert_ne!(t.seed, seed_for(0));
        assert!(t.check_layout().unwrap());
        check_occs(&t, keys, absent, "reseeded");
        // Same input, any order: same file.
        let mut shuffled = occs;
        shuffled.reverse();
        let (mut a, mut b) = (Vec::new(), Vec::new());
        t.save(&mut a).unwrap();
        PoisonTable::build_from_occs(shuffled)
            .unwrap()
            .save(&mut b)
            .unwrap();
        assert_eq!(a, b);
        let t2 = PoisonTable::load(&mut a.as_slice()).unwrap();
        assert_eq!(t2.seed, t.seed);
        check_occs(&t2, keys, absent, "reseeded, reloaded");
    }

    /// A valid `PPOIS02` file whose layout exceeds the probe budget (which
    /// `save` never writes) is re-laid out on load, not rejected or used as is.
    #[test]
    fn test_clustered_file_is_relaid_on_load() {
        let seed = seed_for(0);
        let all = colliding_keys(seed, 3100);
        let (keys, absent) = all.split_at(3000);
        let mut t = PoisonTable::build_from_occs(occs_of(keys)).unwrap();
        let mut entries: Vec<Slot> = t
            .slots
            .iter()
            .copied()
            .filter(|s| *s != Slot::EMPTY)
            .collect();
        entries.sort_unstable_by_key(|s| s.kmer());
        t.slots = PoisonTable::place(&entries, seed, None).unwrap();
        t.seed = seed;
        assert!(
            !t.check_layout().unwrap(),
            "test layout should be clustered"
        );
        let mut bytes = Vec::new();
        t.save(&mut bytes).unwrap();
        let t2 = PoisonTable::load(&mut bytes.as_slice()).unwrap();
        assert_ne!(t2.seed, seed);
        assert!(t2.check_layout().unwrap());
        check_occs(&t2, keys, absent, "re-laid out");
    }

    #[test]
    fn test_corrupt_v2_rejected() {
        let occs = (0..20u64)
            .map(|i| LabeledPoisonOcc {
                canonical_kmer: CanonicalKmer::new(i * 7919),
                unitig_id: i as u32,
                unitig_pos: 1,
            })
            .collect();
        let t = PoisonTable::build_from_occs(occs).unwrap();
        let mut v2 = Vec::new();
        t.save(&mut v2).unwrap();
        // Move a key out of its probe sequence: swap the first occupied slot's
        // key with an empty slot's far away.
        let slot_bytes = &mut v2[8 + 40..];
        let occupied = (0..t.slots.len())
            .find(|&i| t.slots[i] != Slot::EMPTY)
            .unwrap();
        let empty = (0..t.slots.len())
            .rev()
            .find(|&i| t.slots[i] == Slot::EMPTY && i.abs_diff(occupied) > 3)
            .unwrap();
        for j in 0..Slot::BYTES {
            slot_bytes.swap(occupied * Slot::BYTES + j, empty * Slot::BYTES + j);
        }
        assert!(PoisonTable::load(&mut v2.as_slice()).is_err());
        // Truncation.
        let mut v2 = Vec::new();
        t.save(&mut v2).unwrap();
        v2.truncate(v2.len() - 3);
        assert!(PoisonTable::load(&mut v2.as_slice()).is_err());
    }

    /// Real `PPOIS01` files written by the previous implementation must load
    /// into an equivalent table: re-saving as `PPOIS01` reproduces the file
    /// byte for byte. Run with
    /// `POISON_V1_FILES=a.poison,b.poison cargo test --release -- --ignored v1_files`.
    #[test]
    #[ignore]
    fn test_real_v1_files_round_trip() {
        let files = std::env::var("POISON_V1_FILES").expect("set POISON_V1_FILES");
        for path in files.split(',') {
            let bytes = std::fs::read(path).unwrap();
            let t = PoisonTable::load(&mut bytes.as_slice()).unwrap();
            let mut out = Vec::new();
            t.save_v1(&mut out).unwrap();
            assert!(out == bytes, "{path}: PPOIS01 re-save differs");
            eprintln!(
                "{path}: {} k-mers, {} occurrences, round-trips byte for byte",
                t.num_poison_kmers(),
                t.num_poison_occs()
            );
        }
    }

    /// The per-base construction `read_words_of` replaced, kept as its oracle.
    fn read_words_of_per_base(km: CanonicalKmer, k: u32) -> (u64, u64) {
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

    #[test]
    fn test_read_words_of_matches_per_base() {
        let mut x = 0x2545_f491_4f6c_dd1du64;
        for k in 1..=31u32 {
            let mask = (1u64 << (2 * k)) - 1;
            for _ in 0..2000 {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                let km = CanonicalKmer::new(x & mask);
                assert_eq!(
                    ReadKmerFilter::read_words_of(km, k),
                    read_words_of_per_base(km, k),
                    "k={k} km={km:?}"
                );
            }
        }
    }

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
