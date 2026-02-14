//! Poison k-mer scanning for mapping quality filtering.
//!
//! Port of C++ `poison_state_t` with `scan_raw_hits()`. When a poison table
//! is loaded, reads whose mapping intervals contain poison k-mers are
//! discarded to improve specificity.

use crate::index::poison_table::PoisonTable;
use crate::mapping::hit_searcher::SkippingStrategy;
use crate::mapping::hits::FragmentEnd;
use crate::mapping::kmer_value::CanonicalKmer;
use crate::mapping::projected_hits::ProjectedHits;

// ---------------------------------------------------------------------------
// Canonical k-mer helpers
// ---------------------------------------------------------------------------

/// Encode a base as 2-bit value: A=0, C=1, G=2, T=3.
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

    // Full homopolymer: all k bases identical.
    let nuc = kmer & 3;
    let mut homo_mask = 0u64;
    for _ in 0..k {
        homo_mask = (homo_mask << 2) | nuc;
    }
    if kmer == homo_mask {
        return true;
    }

    // XOR of the k-mer with itself shifted one nucleotide left (and nuc
    // filled in at the bottom).  Consecutive identical bases produce zeros.
    // Mask to 2k bits to prevent the topmost nucleotide from overflowing
    // into higher bit positions and corrupting the prefix check.
    let kmask = if k >= 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };
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

/// Iterator over canonical k-mers in a byte sequence.
pub(crate) struct CanonicalKmerIter<'a> {
    seq: &'a [u8],
    k: u32,
    pos: usize,
    fw_kmer: u64,
    rc_kmer: u64,
    mask: u64,
    valid_bases: u32,
}

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
            valid_bases: 0,
        }
    }
}

impl Iterator for CanonicalKmerIter<'_> {
    /// Returns `(canonical_kmer, read_position, is_low_complexity)`.
    type Item = (CanonicalKmer, usize, bool);

    fn next(&mut self) -> Option<Self::Item> {
        while self.pos < self.seq.len() {
            let b = self.seq[self.pos];
            self.pos += 1;

            if let Some(bits) = base_to_bits(b) {
                self.fw_kmer = ((self.fw_kmer << 2) | bits) & self.mask;
                self.rc_kmer = (self.rc_kmer >> 2)
                    | (complement_bits(bits) << (2 * (self.k - 1)));
                self.valid_bases += 1;

                if self.valid_bases >= self.k {
                    let read_pos = self.pos - self.k as usize;
                    let canonical = CanonicalKmer::new(self.fw_kmer.min(self.rc_kmer));
                    let low_complex = is_low_complexity(self.fw_kmer, self.k);
                    return Some((canonical, read_pos, low_complex));
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
    /// Port of C++ `poison_state_t::scan_raw_hits()`.
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

        if raw_hits.is_empty() {
            return false;
        }

        let strict_mode = strat == SkippingStrategy::Strict;
        let first_pos = raw_hits.first().unwrap().0 as usize;

        // Compute terminal position for the last hit.
        let last_phit = &raw_hits.last().unwrap().1;
        let last_uni = last_phit.contig_id();
        let dist_to_contig_end: i64 = if last_phit.hit_fw_on_contig() {
            last_phit.contig_len() as i64 - (last_phit.contig_pos() as i64 + k as i64)
        } else {
            last_phit.contig_pos() as i64
        };
        let terminal_pos = raw_hits.last().unwrap().0 as i64 + dist_to_contig_end;

        let kit = CanonicalKmerIter::new(seq, k);

        // Phase 1: scan k-mers before the first hit.
        for (canonical, read_pos, _low_complex) in kit {
            if read_pos >= first_pos {
                break;
            }
            if ptab.key_exists(canonical) {
                return true;
            }
        }

        // Phase 2: scan intervals between consecutive hits.
        // Re-create iterator for the interval/trailing phases.
        let mut kit = CanonicalKmerIter::new(seq, k);
        // Skip to the first hit position.
        let mut current_kmer: Option<(CanonicalKmer, usize, bool)>;
        loop {
            match kit.next() {
                Some(item) => {
                    if item.1 >= first_pos {
                        current_kmer = Some(item);
                        break;
                    }
                }
                None => return false,
            }
        }

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
                if max_cp < contig_len - k { max_cp - 1 } else { contig_len - k }
            };

            // Scan k-mers in the interval [current_pos, p2)
            loop {
                let item = match current_kmer.take() {
                    Some(i) => i,
                    None => match kit.next() {
                        Some(i) => i,
                        None => break,
                    },
                };
                let (canonical, read_pos, low_complex) = item;

                if read_pos >= p2 {
                    current_kmer = Some((canonical, read_pos, low_complex));
                    break;
                }

                if !low_complex {
                    let was_poisoned = if !right_bound_open {
                        if strict_mode {
                            false
                        } else {
                            ptab.key_occurs_in_unitig_between(canonical, u1, lb, ub)
                        }
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

        // Phase 3: scan remaining k-mers after the last hit interval.
        loop {
            let item = match current_kmer.take() {
                Some(i) => i,
                None => match kit.next() {
                    Some(i) => i,
                    None => break,
                },
            };
            let (canonical, read_pos, _low_complex) = item;
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
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

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
        let kmers: Vec<(CanonicalKmer, usize, bool)> = CanonicalKmerIter::new(seq, 3).collect();
        // "ACGTACGT" has 6 valid 3-mers at positions 0..5
        assert_eq!(kmers.len(), 6);
        assert_eq!(kmers[0].1, 0); // first at position 0
        assert_eq!(kmers[5].1, 5); // last at position 5
    }

    #[test]
    fn test_canonical_kmer_iter_with_n() {
        let seq = b"ACNGT";
        let kmers: Vec<(CanonicalKmer, usize, bool)> = CanonicalKmerIter::new(seq, 3).collect();
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
