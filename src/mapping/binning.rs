//! Genome binning for scATAC-seq mapping.
//!
//! Maps (transcript_id, position) pairs to bin IDs for genomic interval
//! queries. Used when mapping reads to binned reference coordinates.
//!
//! Port of C++ `mapping::util::bin_pos`.

use crate::index::reference_index::ReferenceIndex;

// ---------------------------------------------------------------------------
// BinPos
// ---------------------------------------------------------------------------

/// Genome binning helper.
///
/// Divides the concatenated reference genome into fixed-size bins with
/// optional overlap. Each genomic position maps to one or two bins.
pub struct BinPos {
    /// Cumulative reference lengths (prefix sums).
    cum_ref_lens: Vec<u64>,
    /// Size of each bin in bases.
    _bin_size: u64,
    /// Overlap between adjacent bins in bases.
    bin_overlap: u64,
    /// Effective stride = bin_size - bin_overlap.
    stride: u64,
}

impl BinPos {
    /// Create a new binning scheme from the reference index.
    pub fn new(index: &ReferenceIndex, bin_size: u64, bin_overlap: u64) -> Self {
        let num_refs = index.num_refs();
        let mut cum_ref_lens = Vec::with_capacity(num_refs + 1);
        cum_ref_lens.push(0);
        let mut total = 0u64;
        for i in 0..num_refs {
            total += index.ref_len(i);
            cum_ref_lens.push(total);
        }
        let stride = bin_size.saturating_sub(bin_overlap).max(1);
        Self {
            cum_ref_lens,
            _bin_size: bin_size,
            bin_overlap,
            stride,
        }
    }

    /// Get the bin ID(s) for a (transcript_id, position) pair.
    ///
    /// Returns `(primary_bin, optional_secondary_bin)`. A secondary bin is
    /// returned when the position falls in an overlap region between bins.
    pub fn get_bin_id(&self, tid: u64, pos: u64) -> (u64, Option<u64>) {
        let global_pos = self.cum_ref_lens[tid as usize] + pos;
        let primary = global_pos / self.stride;

        // Check if position falls in overlap region
        let pos_in_bin = global_pos % self.stride;
        let secondary = if self.bin_overlap > 0
            && pos_in_bin < self.bin_overlap
            && primary > 0
        {
            Some(primary - 1)
        } else {
            None
        };

        (primary, secondary)
    }

    /// Total number of bins covering the genome.
    pub fn num_bins(&self) -> u64 {
        let total_len = *self.cum_ref_lens.last().unwrap_or(&0);
        if total_len == 0 {
            return 0;
        }
        total_len.div_ceil(self.stride)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // Helper: create a mock BinPos without needing a real index
    fn make_bin_pos(ref_lens: &[u64], bin_size: u64, bin_overlap: u64) -> BinPos {
        let mut cum = Vec::with_capacity(ref_lens.len() + 1);
        cum.push(0);
        let mut total = 0u64;
        for &len in ref_lens {
            total += len;
            cum.push(total);
        }
        let stride = bin_size.saturating_sub(bin_overlap).max(1);
        BinPos {
            cum_ref_lens: cum,
            _bin_size: bin_size,
            bin_overlap,
            stride,
        }
    }

    #[test]
    fn test_bin_pos_basic() {
        // 2 refs, each 1000bp, bin_size=500, no overlap
        let bp = make_bin_pos(&[1000, 1000], 500, 0);
        assert_eq!(bp.num_bins(), 4); // 2000 / 500 = 4

        let (bin, secondary) = bp.get_bin_id(0, 0);
        assert_eq!(bin, 0);
        assert!(secondary.is_none());

        let (bin, _) = bp.get_bin_id(0, 499);
        assert_eq!(bin, 0);

        let (bin, _) = bp.get_bin_id(0, 500);
        assert_eq!(bin, 1);

        let (bin, _) = bp.get_bin_id(1, 0);
        assert_eq!(bin, 2); // global pos = 1000, bin = 1000/500 = 2
    }

    #[test]
    fn test_bin_pos_overlap_region() {
        // bin_size=1000, bin_overlap=300, stride=700
        let bp = make_bin_pos(&[10000], 1000, 300);
        assert_eq!(bp.stride, 700);

        // Position 0: bin 0, no secondary
        let (bin, secondary) = bp.get_bin_id(0, 0);
        assert_eq!(bin, 0);
        assert!(secondary.is_none()); // bin 0, no previous bin

        // Position 700: bin 1, in overlap (pos_in_bin = 0 < 300), secondary = bin 0
        let (bin, secondary) = bp.get_bin_id(0, 700);
        assert_eq!(bin, 1);
        assert_eq!(secondary, Some(0));

        // Position 800: bin 1, in overlap (pos_in_bin = 100 < 300), secondary = bin 0
        let (bin, secondary) = bp.get_bin_id(0, 800);
        assert_eq!(bin, 1);
        assert_eq!(secondary, Some(0));

        // Position 1000: bin 1, not in overlap (pos_in_bin = 300 >= 300)
        let (bin, secondary) = bp.get_bin_id(0, 1000);
        assert_eq!(bin, 1);
        assert!(secondary.is_none());
    }
}
