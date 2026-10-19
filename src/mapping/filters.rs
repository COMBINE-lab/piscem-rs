//! Poison k-mer scanning for mapping quality filtering.
//!
//! Port of C++ `poison_state_t` with `scan_raw_hits()`. When a poison table
//! is loaded, reads whose mapping intervals contain poison k-mers are
//! discarded to improve specificity.

use crate::index::poison_table::{PoisonTable, QgramSampler, ReadKmerFilter, reverse_groups};
use crate::mapping::hit_searcher::SkippingStrategy;
use crate::mapping::hits::FragmentEnd;
use crate::mapping::kmer_value::CanonicalKmer;
use crate::mapping::projected_hits::ProjectedHits;

// ---------------------------------------------------------------------------
// Canonical k-mer helpers
// ---------------------------------------------------------------------------

/// Encode a base as 2-bit value: A=0, C=1, G=2, T=3.
#[cfg(test)]
#[inline]
fn base_to_bits(b: u8) -> Option<u64> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// Complement of a 2-bit encoded base.
#[cfg(test)]
#[inline]
fn complement_bits(b: u64) -> u64 {
    b ^ 3
}

/// Check if a k-mer (given as 2-bit packed u64) is low-complexity.
///
/// Mirrors C++ `Kmer::is_low_complexity()`: returns true if the k-mer is a
/// full homopolymer, or has a homopolymer run covering at least half of the
/// k-mer at either end (first k/2 bases or last k - k/2 bases in reading
/// order).
///
/// Rust packs k-mers in forward reading order (first char at MSB, last at
/// LSB), so the upper bits correspond to the reading-order prefix and the
/// lower bits to the suffix.  The XOR trick (`kmer ^ ((kmer << 2) | nuc)`)
/// places zero at each position where consecutive bases are identical.  We
/// mask the XOR to 2k bits to avoid overflow at the MSB boundary.
#[inline]
fn is_low_complexity(kmer: u64, k: u32) -> bool {
    if k == 0 {
        return false;
    }

    let kmask = if k >= 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    };

    // Full homopolymer: all k bases identical. `nuc` replicated into each of
    // the k 2-bit groups is `nuc * 0b0101..01` (masked to 2k bits).
    let nuc = kmer & 3;
    let homo_mask = nuc * (0x5555_5555_5555_5555 & kmask);
    if kmer == homo_mask {
        return true;
    }

    // XOR of the k-mer with itself shifted one nucleotide left (and nuc
    // filled in at the bottom).  Consecutive identical bases produce zeros.
    // Mask to 2k bits to prevent the topmost nucleotide from overflowing
    // into higher bit positions and corrupting the prefix check.
    let xor_val = (kmer ^ ((kmer << 2) | nuc)) & kmask;

    // Prefix check: first (k - k/2) reading-order chars all identical.
    // These occupy the upper bit positions in Rust's forward packing.
    let m_prefix = k / 2;
    if (xor_val >> (2 * m_prefix)) == 0 {
        return true;
    }

    // Suffix check: last (k - k/2) reading-order chars all identical.
    // These occupy the lower bit positions; the left shift discards upper bits.
    let m_suffix = k - k / 2;
    if (xor_val << (2 * m_suffix)) == 0 {
        return true;
    }

    false
}

/// Iterator over canonical k-mers in a byte sequence: the original scan's
/// k-mer source, kept for the test oracle.
#[cfg(test)]
pub(crate) struct CanonicalKmerIter<'a> {
    seq: &'a [u8],
    k: u32,
    pos: usize,
    fw_kmer: u64,
    rc_kmer: u64,
    mask: u64,
    /// `2 * (k - 1)`: shift that places a new base at the top of `rc_kmer`.
    rc_shift: u32,
    valid_bases: u32,
}

#[cfg(test)]
impl<'a> CanonicalKmerIter<'a> {
    pub(crate) fn new(seq: &'a [u8], k: u32) -> Self {
        let mask = if k >= 32 {
            u64::MAX
        } else {
            (1u64 << (2 * k)) - 1
        };
        Self {
            seq,
            k,
            pos: 0,
            fw_kmer: 0,
            rc_kmer: 0,
            mask,
            rc_shift: 2 * (k - 1),
            valid_bases: 0,
        }
    }
}

#[cfg(test)]
impl Iterator for CanonicalKmerIter<'_> {
    /// Returns `(canonical_kmer, forward_kmer, read_position)`.
    ///
    /// The forward k-mer is returned (rather than a low-complexity flag) so
    /// callers can run [`is_low_complexity`] only on the rare k-mers whose
    /// verdict matters — those that are actually in the poison table.
    type Item = (CanonicalKmer, u64, usize);

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        while self.pos < self.seq.len() {
            let b = self.seq[self.pos];
            self.pos += 1;

            if let Some(bits) = base_to_bits(b) {
                self.fw_kmer = ((self.fw_kmer << 2) | bits) & self.mask;
                self.rc_kmer = (self.rc_kmer >> 2) | (complement_bits(bits) << self.rc_shift);
                self.valid_bases += 1;

                if self.valid_bases >= self.k {
                    let read_pos = self.pos - self.k as usize;
                    let canonical = CanonicalKmer::new(self.fw_kmer.min(self.rc_kmer));
                    return Some((canonical, self.fw_kmer, read_pos));
                }
            } else {
                // Invalid base — reset
                self.valid_bases = 0;
                self.fw_kmer = 0;
                self.rc_kmer = 0;
            }
        }
        None
    }
}

// ---------------------------------------------------------------------------
// PoisonState
// ---------------------------------------------------------------------------

/// Poison scanning state for a read or read pair.
///
/// Corresponds to C++ `poison_state_t`.
pub struct PoisonState<'a> {
    pub fend: FragmentEnd,
    pub poisoned_left: bool,
    pub poisoned_right: bool,
    pub paired_for_mapping: bool,
    pub ptab: Option<&'a PoisonTable>,
}

impl<'a> PoisonState<'a> {
    /// Create a new poison state with the given (optional) poison table.
    pub fn new(ptab: Option<&'a PoisonTable>) -> Self {
        Self {
            fend: FragmentEnd::Left,
            poisoned_left: false,
            poisoned_right: false,
            paired_for_mapping: false,
            ptab,
        }
    }

    /// Whether poison checking is active.
    #[inline]
    pub fn is_valid(&self) -> bool {
        self.ptab.is_some()
    }

    /// Whether this read (or pair) is poisoned.
    #[inline]
    pub fn is_poisoned(&self) -> bool {
        if self.paired_for_mapping {
            self.poisoned_left || self.poisoned_right
        } else {
            self.poisoned_left
        }
    }

    /// Mark the current fragment end as poisoned.
    #[inline]
    pub fn poison_read(&mut self) {
        match self.fend {
            FragmentEnd::Left => self.poisoned_left = true,
            FragmentEnd::Right => self.poisoned_right = true,
        }
    }

    /// Reset poison state between reads.
    pub fn clear(&mut self) {
        self.poisoned_left = false;
        self.poisoned_right = false;
        self.fend = FragmentEnd::Left;
    }

    /// Set which end of the pair is being processed.
    #[inline]
    pub fn set_fragment_end(&mut self, fend: FragmentEnd) {
        self.fend = fend;
    }

    /// Scan raw hits for poison k-mers.
    ///
    /// Returns `true` if the read is poisoned.
    ///
    /// Port of C++ `poison_state_t::scan_raw_hits()`. That routine walks the
    /// read's k-mers in three phases — before the first hit, within each
    /// interval between consecutive hits, and after the last hit — and
    /// applies a phase-specific table query to every k-mer. Its verdict is
    /// the OR over k-mers, and which phase (and interval) a k-mer falls in
    /// depends only on its read position, so it is computed here only for
    /// k-mers that pass the table's read-side pre-filter; see
    /// [`candidate_is_poisoned`].
    pub fn scan_raw_hits(
        &self,
        seq: &[u8],
        k: u32,
        raw_hits: &[(i32, ProjectedHits<'_>)],
        strat: SkippingStrategy,
    ) -> bool {
        let ptab = match self.ptab {
            Some(pt) => pt,
            None => return false,
        };

        if raw_hits.is_empty() || ptab.empty() {
            return false;
        }

        if seq.len() < k as usize {
            return false;
        }
        let filter = ptab.read_filter(k);
        match &filter.sampler {
            Some(sampler) => scan_sampled(ptab, sampler, seq, k, raw_hits, strat),
            None => scan_rolling(ptab, filter, seq, k, raw_hits, strat),
        }
    }
}

/// Encode a read base as the hit searcher's `ReadKmerIter` does: `(b >> 1) &
/// 3` (A=0 C=1 T=2 G=3). Branchless; non-ACGT bytes map to *some* base, so a
/// window spanning one is rejected later by `candidate_is_poisoned`, which
/// re-reads the window.
#[inline(always)]
fn read_code(b: u8) -> u64 {
    u64::from((b >> 1) & 3)
}

/// Every-position scan: roll the forward k-mer word (base i of the window in
/// bits 2i..2i+2, each new base entering at the top) and probe the read
/// filter at each window. Used for tables too large for [`scan_sampled`].
fn scan_rolling(
    ptab: &PoisonTable,
    filter: &ReadKmerFilter,
    seq: &[u8],
    k: u32,
    raw_hits: &[(i32, ProjectedHits<'_>)],
    strat: SkippingStrategy,
) -> bool {
    let ku = k as usize;
    let top_shift = 2 * (k - 1);
    let roll = |fw: u64, b: u8| (fw >> 2) | (read_code(b) << top_shift);
    let mut fw = seq[..ku - 1].iter().fold(0u64, |fw, &b| roll(fw, b));
    for (pos, &b) in seq[ku - 1..].iter().enumerate() {
        fw = roll(fw, b);
        if filter.may_contain(fw)
            && candidate_is_poisoned(ptab, &seq[pos..pos + ku], fw, pos, k, raw_hits, strat)
        {
            return true;
        }
    }
    false
}

/// Read-encoded code of the `q <= 16` bases at `seq[at..at + q]`, base i in
/// bits 2i..2i+2 (bits above 2q may be set; callers mask).
#[inline(always)]
fn qgram_code(seq: &[u8], at: usize, q: usize) -> u64 {
    // Pack the 2-bit codes of 8 bytes (byte j's code in bits 8j..8j+2) into
    // 16 bits, base j at bits 2j..2j+2.
    #[inline(always)]
    fn pack8(bytes: u64) -> u64 {
        let c = (bytes >> 1) & 0x0303_0303_0303_0303;
        let c = (c | (c >> 6)) & 0x000F_000F_000F_000F;
        let c = (c | (c >> 12)) & 0x0000_00FF_0000_00FF;
        (c | (c >> 24)) & 0xFFFF
    }
    let load8 = |i: usize| u64::from_le_bytes(seq[i..i + 8].try_into().unwrap());
    if q <= 8 && at + 8 <= seq.len() {
        pack8(load8(at))
    } else if at + 16 <= seq.len() {
        pack8(load8(at)) | (pack8(load8(at + 8)) << 16)
    } else {
        seq[at..at + q]
            .iter()
            .enumerate()
            .fold(0u64, |c, (i, &b)| c | (read_code(b) << (2 * i)))
    }
}

/// Sampled scan for small tables: test the q-gram at every `stride`-th read
/// position against the sampler's bitmap, and confirm only the windows a hit
/// implicates. Every window that is a poison k-mer contains a sampled q-gram
/// that the bitmap knows (see [`QgramSampler`]), so the verdict matches
/// [`scan_rolling`]'s.
fn scan_sampled(
    ptab: &PoisonTable,
    sampler: &QgramSampler,
    seq: &[u8],
    k: u32,
    raw_hits: &[(i32, ProjectedHits<'_>)],
    strat: SkippingStrategy,
) -> bool {
    let ku = k as usize;
    let q = sampler.q as usize;
    let n_kmers = seq.len() + 1 - ku;
    let qmask = (1u64 << (2 * q)) - 1;
    let mut sample = 0usize;
    while sample + q <= seq.len() {
        let code = (qgram_code(seq, sample, q) & qmask) as u32;
        if sampler.has(code) {
            for (offset, pattern) in sampler.patterns_with(code) {
                let Some(pos) = sample.checked_sub(offset as usize) else {
                    continue;
                };
                if pos >= n_kmers {
                    continue;
                }
                let window = &seq[pos..pos + ku];
                let word = window
                    .iter()
                    .enumerate()
                    .fold(0u64, |w, (i, &b)| w | (read_code(b) << (2 * i)));
                if word == pattern
                    && candidate_is_poisoned(ptab, window, word, pos, k, raw_hits, strat)
                {
                    return true;
                }
            }
        }
        sample += sampler.stride;
    }
    false
}

/// The canonical k-mer and forward k-mer (both ACGT, first base most
/// significant — the poison table's key layout) of the window `kmer`,
/// given its read-encoded forward word `word`; `None` if the window holds a
/// non-ACGT byte (which the read encoding silently maps to some base).
#[inline]
fn canonical_from_read_word(kmer: &[u8], word: u64, k: u32) -> Option<(CanonicalKmer, u64)> {
    if !kmer
        .iter()
        .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't'))
    {
        return None;
    }
    // Read code -> ACGT code per 2-bit group (swap G and T); base i stays in
    // bits 2i..2i+2. Reversing the group order puts the first base on top;
    // complementing (`^ 0b11` per group) without reversing is the reverse
    // complement in that same first-base-on-top layout.
    let c = word ^ ((word >> 1) & 0x5555_5555_5555_5555);
    let fw = reverse_groups(c, k);
    let rc = c ^ ((1u64 << (2 * k)) - 1);
    Some((CanonicalKmer::new(fw.min(rc)), fw))
}

/// The exact poison verdict for the k-mer `kmer` at read position `pos`.
///
/// Applies the query `scan_raw_hits()`'s phase for `pos` would apply:
///
/// * before the first hit: the k-mer must not be poison anywhere;
/// * in an interval between consecutive hits — the first pair `(h_j, h_j+1)`
///   with `pos < pos(h_j+1)`, i.e. the pair whose sequential scan would
///   consume this k-mer — poison is checked on `h_j`'s unitig between the
///   two hits' unitig positions (or anywhere, if `h_j+1` came from an open
///   search; never, under strict skipping with a closed interval), and
///   low-complexity k-mers are exempt;
/// * after the last hit: poison is checked on the last hit's unitig while the
///   k-mer could still lie on it, and anywhere after that.
///
/// Only reached for pre-filter hits, so kept out of the scan loop.
#[cold]
#[inline(never)]
fn candidate_is_poisoned(
    ptab: &PoisonTable,
    kmer: &[u8],
    word: u64,
    pos: usize,
    k: u32,
    raw_hits: &[(i32, ProjectedHits<'_>)],
    strat: SkippingStrategy,
) -> bool {
    // `None` iff the window holds a non-ACGT byte: no such k-mer exists.
    let Some((canonical, fw_kmer)) = canonical_from_read_word(kmer, word, k) else {
        return false;
    };

    // Phase 1: before the first hit.
    let first_pos = raw_hits[0].0 as usize;
    if pos < first_pos {
        return ptab.key_exists(canonical);
    }

    // Phase 2: the interval whose sequential scan would consume `pos`.
    if let Some(pair) = raw_hits.windows(2).find(|p| pos < p[1].0 as usize) {
        let (start_hit, end_hit) = (&pair[0].1, &pair[1].1);
        let right_bound_open = end_hit.resulted_from_open_search();
        if !right_bound_open && strat == SkippingStrategy::Strict {
            return false;
        }
        let in_table = if right_bound_open {
            ptab.key_exists(canonical)
        } else {
            let u1 = start_hit.contig_id();
            let cp1 = start_hit.contig_pos();
            let cp2 = end_hit.contig_pos();
            let contig_len = start_hit.contig_len();
            let lb = {
                let min_cp = cp1.min(cp2);
                if min_cp > 0 { min_cp + 1 } else { 0 }
            };
            let ub = {
                let max_cp = cp1.max(cp2);
                if max_cp < contig_len - k {
                    max_cp - 1
                } else {
                    contig_len - k
                }
            };
            ptab.key_occurs_in_unitig_between(canonical, u1, lb, ub)
        };
        return in_table && !is_low_complexity(fw_kmer, k);
    }

    // Phase 3: after the last hit.
    let (last_read_pos, last_phit) = raw_hits.last().unwrap();
    let dist_to_contig_end: i64 = if last_phit.hit_fw_on_contig() {
        last_phit.contig_len() as i64 - (last_phit.contig_pos() as i64 + k as i64)
    } else {
        last_phit.contig_pos() as i64
    };
    let terminal_pos = *last_read_pos as i64 + dist_to_contig_end;
    if (pos as i64) < terminal_pos {
        ptab.key_occurs_in_unitig(canonical, last_phit.contig_id())
    } else {
        ptab.key_exists(canonical)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::index::contig_table::ContigTableBuilder;
    use crate::index::poison_table::LabeledPoisonOcc;
    use crate::mapping::hit_searcher::ReadKmerIter;

    /// The original three-phase scan (a direct port of C++
    /// `poison_state_t::scan_raw_hits()`), kept as the oracle for
    /// `scan_raw_hits`. Visits every k-mer and computes low-complexity
    /// eagerly, exactly as the port did.
    fn reference_scan(
        ptab: &PoisonTable,
        seq: &[u8],
        k: u32,
        raw_hits: &[(i32, ProjectedHits<'_>)],
        strat: SkippingStrategy,
    ) -> bool {
        if raw_hits.is_empty() {
            return false;
        }
        let strict_mode = strat == SkippingStrategy::Strict;
        let first_pos = raw_hits.first().unwrap().0 as usize;
        let last_phit = &raw_hits.last().unwrap().1;
        let last_uni = last_phit.contig_id();
        let dist_to_contig_end: i64 = if last_phit.hit_fw_on_contig() {
            last_phit.contig_len() as i64 - (last_phit.contig_pos() as i64 + k as i64)
        } else {
            last_phit.contig_pos() as i64
        };
        let terminal_pos = raw_hits.last().unwrap().0 as i64 + dist_to_contig_end;

        let items = |seq| {
            CanonicalKmerIter::new(seq, k).map(move |(c, fw, p)| (c, p, is_low_complexity(fw, k)))
        };
        for (canonical, read_pos, _) in items(seq) {
            if read_pos >= first_pos {
                break;
            }
            if ptab.key_exists(canonical) {
                return true;
            }
        }
        let mut kit = items(seq);
        let mut current_kmer = loop {
            match kit.next() {
                Some(item) if item.1 >= first_pos => break Some(item),
                Some(_) => {}
                None => return false,
            }
        };
        let mut start_idx = 0usize;
        let mut end_idx = 1usize;
        while end_idx < raw_hits.len() {
            let u1 = raw_hits[start_idx].1.contig_id();
            let cp1 = raw_hits[start_idx].1.contig_pos();
            let cp2 = raw_hits[end_idx].1.contig_pos();
            let p2 = raw_hits[end_idx].0 as usize;
            let right_bound_open = raw_hits[end_idx].1.resulted_from_open_search();
            let contig_len = raw_hits[start_idx].1.contig_len();
            let lb = {
                let min_cp = cp1.min(cp2);
                if min_cp > 0 { min_cp + 1 } else { 0 }
            };
            let ub = {
                let max_cp = cp1.max(cp2);
                if max_cp < contig_len - k {
                    max_cp - 1
                } else {
                    contig_len - k
                }
            };
            while let Some(item) = current_kmer.take().or_else(|| kit.next()) {
                let (canonical, read_pos, low_complex) = item;
                if read_pos >= p2 {
                    current_kmer = Some(item);
                    break;
                }
                if !low_complex {
                    let was_poisoned = if !right_bound_open {
                        !strict_mode && ptab.key_occurs_in_unitig_between(canonical, u1, lb, ub)
                    } else {
                        ptab.key_exists(canonical)
                    };
                    if was_poisoned {
                        return true;
                    }
                }
            }
            start_idx += 1;
            end_idx += 1;
        }
        while let Some((canonical, read_pos, _)) = current_kmer.take().or_else(|| kit.next()) {
            let was_poisoned = if (read_pos as i64) < terminal_pos {
                ptab.key_occurs_in_unitig(canonical, last_uni)
            } else {
                ptab.key_exists(canonical)
            };
            if was_poisoned {
                return true;
            }
        }
        false
    }

    /// Deterministic xorshift64* for the randomized tests.
    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 >> 12;
            self.0 ^= self.0 << 25;
            self.0 ^= self.0 >> 27;
            self.0.wrapping_mul(0x2545_F491_4F6C_DD1D)
        }
        fn below(&mut self, n: u64) -> u64 {
            self.next() % n
        }
    }

    fn canonical_of(kmer: &[u8], k: u32) -> CanonicalKmer {
        CanonicalKmerIter::new(kmer, k).next().unwrap().0
    }

    fn revcomp(s: &[u8]) -> Vec<u8> {
        s.iter()
            .rev()
            .map(|&b| match b.to_ascii_uppercase() {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                x => x,
            })
            .collect()
    }

    /// The read-side filter must hold each poison k-mer under the exact words
    /// the hit searcher's `ReadKmerIter` produces for it, in both
    /// orientations — otherwise the scan would miss poison k-mers.
    #[test]
    fn test_read_filter_matches_read_kmer_iter_encoding() {
        let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
        for k in [5u32, 19, 23, 31] {
            let kmers: Vec<Vec<u8>> = (0..2000)
                .map(|_| (0..k).map(|_| b"ACGTacgt"[rng.below(8) as usize]).collect())
                .collect();
            let occs = kmers
                .iter()
                .map(|s| LabeledPoisonOcc {
                    canonical_kmer: canonical_of(s, k),
                    unitig_id: 0,
                    unitig_pos: 0,
                })
                .collect();
            let pt = PoisonTable::build_from_occs(occs).unwrap();
            let filter = pt.read_filter(k);
            for s in &kmers {
                for oriented in [s.clone(), revcomp(s)] {
                    let word = ReadKmerIter::new(&oriented, k as usize).fw_word();
                    assert!(
                        filter.may_contain(word),
                        "k={k}: filter misses {:?}",
                        oriented
                    );
                }
            }
        }
    }

    #[test]
    fn test_qgram_code_matches_per_base() {
        let mut rng = Rng(0x51_7cc1_b727_220a);
        let seq: Vec<u8> = (0..97)
            .map(|_| b"ACGTacgtN"[rng.below(9) as usize])
            .collect();
        for q in 1..=16usize {
            let mask = (1u64 << (2 * q)) - 1;
            for at in 0..=seq.len() - q {
                let expected = seq[at..at + q]
                    .iter()
                    .enumerate()
                    .fold(0u64, |c, (i, &b)| c | (read_code(b) << (2 * i)));
                assert_eq!(qgram_code(&seq, at, q) & mask, expected, "q={q} at={at}");
            }
        }
    }

    #[test]
    fn test_canonical_from_read_word_matches_iter() {
        let mut rng = Rng(0x6c62_272e_07bb_0142);
        for k in [1u32, 5, 19, 23, 31] {
            for _ in 0..3000 {
                let window: Vec<u8> = (0..k)
                    .map(|_| match rng.below(40) {
                        0 => b'N',
                        _ => b"ACGTacgt"[rng.below(8) as usize],
                    })
                    .collect();
                let word = window
                    .iter()
                    .enumerate()
                    .fold(0u64, |w, (i, &b)| w | (read_code(b) << (2 * i)));
                let expected = CanonicalKmerIter::new(&window, k)
                    .next()
                    .map(|(c, fw, _)| (c, fw));
                assert_eq!(
                    canonical_from_read_word(&window, word, k),
                    expected,
                    "k={k} window={}",
                    String::from_utf8_lossy(&window)
                );
            }
        }
    }

    /// `scan_raw_hits` must give the original scan's verdict on every input.
    #[test]
    fn test_scan_raw_hits_matches_reference() {
        let mut builder = ContigTableBuilder::new(1, 10_000, 1);
        builder.add_occurrence(0, 0, 100, true);
        let table = builder.build();

        let mut rng = Rng(0x2545_f491_4f6c_dd1d);
        let mut n_poisoned = 0usize;
        let n_trials = 20_000;
        for trial in 0..n_trials {
            let k = [5u32, 7, 11, 23][trial % 4];
            let len = k as usize + rng.below(60) as usize;
            let mut read: Vec<u8> = (0..len)
                .map(|_| match rng.below(100) {
                    0 => b'N',
                    1 => b'a',
                    _ => b"ACGT"[rng.below(4) as usize],
                })
                .collect();
            // Occasionally plant a homopolymer run to exercise the
            // low-complexity exemption.
            if rng.below(4) == 0 {
                let start = rng.below(len as u64) as usize;
                let end = (start + k as usize).min(len);
                let base = b"ACGT"[rng.below(4) as usize];
                read[start..end].fill(base);
            }

            // Poison k-mers: some read windows (so hits happen), some random.
            let n_unitigs = 3u32;
            let mut occs = Vec::new();
            for _ in 0..1 + rng.below(4) {
                let canonical = if rng.below(3) > 0 && len >= k as usize {
                    let p = rng.below((len - k as usize + 1) as u64) as usize;
                    match CanonicalKmerIter::new(&read[p..p + k as usize], k).next() {
                        Some(item) => item.0,
                        None => continue,
                    }
                } else {
                    CanonicalKmer::new(rng.next() & ((1u64 << (2 * k)) - 1))
                };
                occs.push(LabeledPoisonOcc {
                    canonical_kmer: canonical,
                    unitig_id: rng.below(n_unitigs as u64) as u32,
                    unitig_pos: rng.below(40) as u32,
                });
            }
            let pt = PoisonTable::build_from_occs(occs).unwrap();

            // Raw hits at non-decreasing read positions, as the hit searcher
            // emits them.
            let n_hits = 1 + rng.below(5) as usize;
            let mut positions: Vec<i32> = (0..n_hits)
                .map(|_| rng.below((len - k as usize + 1) as u64) as i32)
                .collect();
            positions.sort_unstable();
            let raw_hits: Vec<(i32, ProjectedHits<'_>)> = positions
                .into_iter()
                .map(|p| {
                    let contig_len = k + 1 + rng.below(60) as u32;
                    let contig_pos = 1 + rng.below((contig_len - k) as u64) as u32;
                    let mut ph = ProjectedHits::new(
                        rng.below(n_unitigs as u64) as u32,
                        contig_pos,
                        rng.below(2) == 0,
                        contig_len,
                        0,
                        k,
                        table.contig_entries(0),
                    );
                    ph.set_resulted_from_open_search(rng.below(3) == 0);
                    (p, ph)
                })
                .collect();

            let ps = PoisonState::new(Some(&pt));
            // Each scan path on its own, the sampled one at a random q-gram
            // length (q = 1..=min(k, 8): every stride down to 1 for small k).
            let filter = pt.read_filter(k);
            let patterns: Vec<u64> = pt
                .keys()
                .flat_map(|km| {
                    let (fw, rc) = ReadKmerFilter::read_words_of(km, k);
                    [fw, rc]
                })
                .collect();
            let sampler =
                QgramSampler::build(&patterns, k, 1 + rng.below(u64::from(k.min(8))) as u32);
            for strat in [SkippingStrategy::Permissive, SkippingStrategy::Strict] {
                let expected = reference_scan(&pt, &read, k, &raw_hits, strat);
                let ctx = || {
                    format!(
                        "trial {trial} k={k} q={} strat={strat:?} read={}",
                        sampler.q,
                        String::from_utf8_lossy(&read)
                    )
                };
                let got = ps.scan_raw_hits(&read, k, &raw_hits, strat);
                assert_eq!(got, expected, "dispatch: {}", ctx());
                let got = scan_rolling(&pt, filter, &read, k, &raw_hits, strat);
                assert_eq!(got, expected, "rolling: {}", ctx());
                let got = scan_sampled(&pt, &sampler, &read, k, &raw_hits, strat);
                assert_eq!(got, expected, "sampled: {}", ctx());
                n_poisoned += usize::from(expected);
            }
        }
        // Both verdicts must be well represented for the test to mean much.
        assert!(n_poisoned > n_trials / 10 && n_poisoned < 2 * n_trials - n_trials / 10);
    }

    #[test]
    fn test_poison_clear() {
        let mut ps = PoisonState::new(None);
        ps.poisoned_left = true;
        ps.poisoned_right = true;
        ps.fend = FragmentEnd::Right;
        ps.clear();
        assert!(!ps.poisoned_left);
        assert!(!ps.poisoned_right);
        assert_eq!(ps.fend, FragmentEnd::Left);
    }

    #[test]
    fn test_no_poison_table_returns_false() {
        let ps = PoisonState::new(None);
        assert!(!ps.is_valid());
        let result = ps.scan_raw_hits(b"ACGTACGT", 5, &[], SkippingStrategy::Permissive);
        assert!(!result);
    }

    #[test]
    fn test_empty_hits_not_poisoned() {
        use crate::index::poison_table::{LabeledPoisonOcc, PoisonTable};
        let occs = vec![LabeledPoisonOcc {
            canonical_kmer: CanonicalKmer::new(42),
            unitig_id: 0,
            unitig_pos: 5,
        }];
        let pt = PoisonTable::build_from_occs(occs).unwrap();
        let ps = PoisonState::new(Some(&pt));
        assert!(ps.is_valid());
        let result = ps.scan_raw_hits(b"ACGTACGT", 5, &[], SkippingStrategy::Permissive);
        assert!(!result);
    }

    #[test]
    fn test_is_poisoned_paired() {
        let mut ps = PoisonState::new(None);
        ps.paired_for_mapping = true;
        assert!(!ps.is_poisoned());

        ps.poisoned_left = true;
        assert!(ps.is_poisoned());

        ps.poisoned_left = false;
        ps.poisoned_right = true;
        assert!(ps.is_poisoned());
    }

    #[test]
    fn test_is_poisoned_single() {
        let mut ps = PoisonState::new(None);
        ps.paired_for_mapping = false;
        assert!(!ps.is_poisoned());

        ps.poisoned_right = true; // Shouldn't matter for single-end
        assert!(!ps.is_poisoned());

        ps.poisoned_left = true;
        assert!(ps.is_poisoned());
    }

    #[test]
    fn test_canonical_kmer_iter() {
        let seq = b"ACGTACGT";
        let kmers: Vec<(CanonicalKmer, u64, usize)> = CanonicalKmerIter::new(seq, 3).collect();
        // "ACGTACGT" has 6 valid 3-mers at positions 0..5
        assert_eq!(kmers.len(), 6);
        assert_eq!(kmers[0].2, 0); // first at position 0
        assert_eq!(kmers[5].2, 5); // last at position 5
    }

    #[test]
    fn test_canonical_kmer_iter_with_n() {
        let seq = b"ACNGT";
        let kmers: Vec<(CanonicalKmer, u64, usize)> = CanonicalKmerIter::new(seq, 3).collect();
        // "ACN" invalid, "CNG" invalid, "NGT" invalid — no valid 3-mers
        // Actually: AC is valid, then N resets. Then G,T only 2 bases. So 0 k-mers.
        assert_eq!(kmers.len(), 0);
    }

    #[test]
    fn test_is_low_complexity_homopolymer() {
        // Full homopolymer: TTTTT (k=5) → low complexity
        let k = 5u32;
        // Rust forward order: T at MSB → T=3, packed = 11_11_11_11_11
        let poly_t: u64 = (0..k).fold(0u64, |acc, _| (acc << 2) | 3);
        assert!(is_low_complexity(poly_t, k));
    }

    #[test]
    fn test_is_low_complexity_prefix_run() {
        // K-mer with a long run of T at the START followed by different bases.
        // This should be detected as low-complexity (matching C++ behavior).
        // Build "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTGC" (29 T's + GC, k=31)
        let k = 31u32;
        let seq = b"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTGC";
        let mut fw = 0u64;
        let mask = (1u64 << (2 * k)) - 1;
        for &b in seq.iter() {
            let bits = base_to_bits(b).unwrap();
            fw = ((fw << 2) | bits) & mask;
        }
        // First 29 chars are T, k/2 = 15, so first 16 chars are all T.
        // C++ suffix check fires → low complexity.
        assert!(is_low_complexity(fw, k));
    }

    #[test]
    fn test_is_low_complexity_fam76a_kmer() {
        // The specific k-mer that caused the parity difference:
        // TTTTTTTTTTTTTTTTTTTTTTGTGGTGGTG (22 T's + GTGGTGGTG, k=31)
        // First 16 chars are all T → C++ marks this as low-complexity.
        let k = 31u32;
        let seq = b"TTTTTTTTTTTTTTTTTTTTTTGTGGTGGTG";
        let mut fw = 0u64;
        let mask = (1u64 << (2 * k)) - 1;
        for &b in seq.iter() {
            let bits = base_to_bits(b).unwrap();
            fw = ((fw << 2) | bits) & mask;
        }
        assert!(is_low_complexity(fw, k));
    }

    #[test]
    fn test_is_low_complexity_suffix_run() {
        // K-mer with different bases at start, then a long run at the END.
        // Build "GCTTTTTTTTTTTTTTTTTTTTTTTTTTTTT" (GC + 29 T's, k=31)
        // The end has 29 T's, last k-k/2=16 chars are all T → should fire.
        let k = 31u32;
        let seq = b"GCTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
        assert_eq!(seq.len(), k as usize);
        let mut fw = 0u64;
        let mask = (1u64 << (2 * k)) - 1;
        for &b in seq.iter() {
            let bits = base_to_bits(b).unwrap();
            fw = ((fw << 2) | bits) & mask;
        }
        assert!(is_low_complexity(fw, k));
    }
}
