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

/// Genome binning helper matching C++ `bin_pos`.
///
/// Divides each reference into fixed-size bins. Bin IDs are cumulative
/// across references: reference 0 gets bins 0..N0-1, reference 1 gets
/// bins N0..N0+N1-1, etc.
pub struct BinPos {
    /// cum_bin_ids[i] = cumulative number of bins for refs 0..i.
    /// Length = num_refs + 1.
    cum_bin_ids: Vec<u64>,
    /// Size of each bin in bases.
    bin_size: u64,
    /// Overlap between adjacent bins in bases.
    overlap: u64,
    /// Hit threshold fraction (e.g., 0.7).
    thr: f32,
}

/// Sentinel for an invalid bin ID.
pub const INVALID_BIN_ID: u64 = u64::MAX;

impl BinPos {
    /// Create a new binning scheme from the reference index.
    ///
    /// Matches C++ `bin_pos::compute_cum_rank()`.
    pub fn new(index: &ReferenceIndex, bin_size: u64, overlap: u64, thr: f32) -> Self {
        let num_refs = index.num_refs();
        let mut cum_bin_ids = Vec::with_capacity(num_refs + 1);
        cum_bin_ids.push(0);
        for i in 0..num_refs {
            let rlen = index.ref_len(i);
            let mut bins_in_ref = rlen / bin_size;
            if bins_in_ref * bin_size < rlen {
                bins_in_ref += 1;
            }
            cum_bin_ids.push(cum_bin_ids[i] + bins_in_ref);
        }
        Self {
            cum_bin_ids,
            bin_size,
            overlap,
            thr,
        }
    }

    /// Get the bin ID(s) for a (tid, pos) pair.
    ///
    /// Returns `(primary_bin, secondary_bin)` where secondary is
    /// `INVALID_BIN_ID` if the position is not in an overlap region.
    ///
    /// Matches C++ `bin_pos::get_bin_id()`.
    #[inline]
    pub fn get_bin_id(&self, tid: u32, pos: u64) -> (u64, u64) {
        let first_bin_id = self.cum_bin_ids[tid as usize];
        let first_bin_id_in_next = self.cum_bin_ids[tid as usize + 1];
        debug_assert!(self.bin_size > self.overlap);

        let rel_bin = pos / self.bin_size;
        let bin1 = first_bin_id + rel_bin;

        let bin2_on_same_ref = (bin1 + 1) < first_bin_id_in_next;
        let bin2_rel_start_pos = (rel_bin + 1) * self.bin_size;
        let bin2 = if bin2_on_same_ref && pos > bin2_rel_start_pos.saturating_sub(self.overlap) {
            bin1 + 1
        } else {
            INVALID_BIN_ID
        };

        (bin1, bin2)
    }

    /// Whether a bin ID is valid (not the sentinel).
    #[inline]
    pub fn is_valid(bin_id: u64) -> bool {
        bin_id != INVALID_BIN_ID
    }

    /// Hit threshold fraction.
    #[inline]
    pub fn thr(&self) -> f32 {
        self.thr
    }

    /// Bin size in bases.
    #[inline]
    pub fn bin_size(&self) -> u64 {
        self.bin_size
    }

    /// Overlap in bases.
    #[inline]
    pub fn overlap(&self) -> u64 {
        self.overlap
    }

    /// Total number of bins across all references.
    pub fn num_bins(&self) -> u64 {
        *self.cum_bin_ids.last().unwrap_or(&0)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // Helper: create a mock BinPos without needing a real index
    fn make_bin_pos(ref_lens: &[u64], bin_size: u64, overlap: u64, thr: f32) -> BinPos {
        let mut cum = Vec::with_capacity(ref_lens.len() + 1);
        cum.push(0u64);
        for &rlen in ref_lens {
            let mut bins_in_ref = rlen / bin_size;
            if bins_in_ref * bin_size < rlen {
                bins_in_ref += 1;
            }
            cum.push(cum.last().unwrap() + bins_in_ref);
        }
        BinPos {
            cum_bin_ids: cum,
            bin_size,
            overlap,
            thr,
        }
    }

    #[test]
    fn test_bin_pos_basic() {
        // 2 refs: 1000bp each, bin_size=500, no overlap
        let bp = make_bin_pos(&[1000, 1000], 500, 0, 0.7);
        assert_eq!(bp.num_bins(), 4); // 2 + 2

        let (bin, sec) = bp.get_bin_id(0, 0);
        assert_eq!(bin, 0);
        assert_eq!(sec, INVALID_BIN_ID);

        let (bin, _) = bp.get_bin_id(0, 499);
        assert_eq!(bin, 0);

        let (bin, _) = bp.get_bin_id(0, 500);
        assert_eq!(bin, 1);

        let (bin, _) = bp.get_bin_id(1, 0);
        assert_eq!(bin, 2);

        let (bin, _) = bp.get_bin_id(1, 500);
        assert_eq!(bin, 3);
    }

    #[test]
    fn test_bin_pos_overlap_region() {
        // bin_size=1000, overlap=300
        // Ref 0: 5000bp → 5 bins (0..4)
        let bp = make_bin_pos(&[5000], 1000, 300, 0.7);
        assert_eq!(bp.num_bins(), 5);

        // Position 0: bin 0, no secondary
        let (bin, sec) = bp.get_bin_id(0, 0);
        assert_eq!(bin, 0);
        assert_eq!(sec, INVALID_BIN_ID);

        // Position 699: bin 0, no secondary (699 <= 1000-300=700, no: 699 > 700? No.)
        let (bin, sec) = bp.get_bin_id(0, 699);
        assert_eq!(bin, 0);
        assert_eq!(sec, INVALID_BIN_ID);

        // Position 701: bin 0, secondary = bin 1 (701 > 1000-300=700)
        let (bin, sec) = bp.get_bin_id(0, 701);
        assert_eq!(bin, 0);
        assert_eq!(sec, 1);

        // Position 999: bin 0, secondary = bin 1
        let (bin, sec) = bp.get_bin_id(0, 999);
        assert_eq!(bin, 0);
        assert_eq!(sec, 1);

        // Position 1000: bin 1, no secondary (1000 <= 2000-300=1700, not > 1700)
        let (bin, sec) = bp.get_bin_id(0, 1000);
        assert_eq!(bin, 1);
        assert_eq!(sec, INVALID_BIN_ID);
    }

    #[test]
    fn test_bin_pos_multi_ref() {
        // Ref 0: 2500bp → ceil(2500/1000) = 3 bins (0,1,2)
        // Ref 1: 1500bp → ceil(1500/1000) = 2 bins (3,4)
        let bp = make_bin_pos(&[2500, 1500], 1000, 300, 0.7);
        assert_eq!(bp.num_bins(), 5);

        // Ref 1, pos 0 → bin 3
        let (bin, _) = bp.get_bin_id(1, 0);
        assert_eq!(bin, 3);

        // Ref 1, pos 1000 → bin 4
        let (bin, _) = bp.get_bin_id(1, 1000);
        assert_eq!(bin, 4);
    }

    #[test]
    fn test_bin_pos_last_ref_boundary() {
        // Ref with exactly bin_size length: 1 bin
        let bp = make_bin_pos(&[1000], 1000, 300, 0.7);
        assert_eq!(bp.num_bins(), 1);

        // Position 0: bin 0, no secondary (no next bin)
        let (bin, sec) = bp.get_bin_id(0, 0);
        assert_eq!(bin, 0);
        assert_eq!(sec, INVALID_BIN_ID);

        // Position 900: bin 0, no secondary (bin1+1=1 not < 1)
        let (bin, sec) = bp.get_bin_id(0, 900);
        assert_eq!(bin, 0);
        assert_eq!(sec, INVALID_BIN_ID);
    }
}
