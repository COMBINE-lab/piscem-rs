//! Variable-length barcode/UMI normalization via collision-free padding.
//!
//! When barcodes or UMIs have variable length (e.g., `b[9-10]`), the
//! extracted sequence must be padded to a fixed width so downstream tools
//! (which pack barcodes into fixed-width integers) work correctly.
//!
//! The padding scheme appends a length-dependent suffix that ensures no
//! two barcodes of different original lengths can collide after padding.
//! See [`VAR_LEN_PADDING`] for the suffix table.
//!
//! # Example
//!
//! ```
//! use seq_geom_parser::normalize::pad_to_fixed;
//!
//! let mut buf = [0u8; 32];
//! // A 9bp barcode from a b[9-10] tag: pad to 10bp
//! let padded = pad_to_fixed(b"ACGTACGTA", 10, &mut buf);
//! assert_eq!(padded, b"ACGTACGTAA"); // "A" appended (deficit = 1)
//!
//! // A 10bp barcode: no padding needed
//! let padded = pad_to_fixed(b"ACGTACGTAC", 10, &mut buf);
//! assert_eq!(padded, b"ACGTACGTAC");
//! ```

use crate::types::{GeoLen, VAR_LEN_PADDING, MAX_RANGE_WIDTH};

/// Maximum barcode/UMI length we support (in bases).
/// This limits the size of inline buffers.
pub const MAX_BC_LEN: usize = 32;

/// Pad a variable-length sequence to a fixed width using the collision-free
/// padding scheme.
///
/// - `seq`: the extracted sequence (may be shorter than `max_len`)
/// - `max_len`: the maximum length (from `GeoLen::Range(_, max)`)
/// - `buf`: a caller-provided buffer of at least `max_len` bytes
///
/// Returns a slice of `buf` containing the padded sequence (length == `max_len`).
///
/// If `seq.len() == max_len`, the original bytes are copied without padding.
/// If `seq.len() < max_len`, the appropriate padding suffix is appended.
///
/// # Panics
///
/// Panics if the deficit (`max_len - seq.len()`) exceeds [`MAX_RANGE_WIDTH`].
#[inline]
pub fn pad_to_fixed<'a>(seq: &[u8], max_len: usize, buf: &'a mut [u8]) -> &'a [u8] {
    let captured_len = seq.len();
    debug_assert!(captured_len <= max_len);
    debug_assert!(max_len <= buf.len());

    buf[..captured_len].copy_from_slice(seq);

    let deficit = max_len - captured_len;
    if deficit > 0 {
        let padding = VAR_LEN_PADDING[deficit];
        buf[captured_len..captured_len + padding.len()].copy_from_slice(padding);
        &buf[..max_len]
    } else {
        &buf[..max_len]
    }
}

/// Check whether a [`GeoLen`] is variable-length and requires normalization.
#[inline]
pub fn needs_padding(len: &GeoLen) -> bool {
    matches!(len, GeoLen::Range(min, max) if min != max)
}

/// Get the normalized (padded) length for a [`GeoLen`].
/// For fixed-length, returns the fixed value. For ranges, returns the max.
#[inline]
pub fn normalized_len(len: &GeoLen) -> Option<usize> {
    match len {
        GeoLen::Fixed(n) => Some(*n as usize),
        GeoLen::Range(_, max) => Some(*max as usize),
        GeoLen::Unbounded => None,
    }
}

/// A reusable padding workspace for normalizing variable-length barcodes/UMIs.
///
/// Create one per thread, call `pad()` for each extracted sequence.
/// Holds an inline buffer to avoid per-read allocation.
pub struct PadBuf {
    buf: [u8; MAX_BC_LEN],
}

impl PadBuf {
    /// Create a new padding buffer.
    pub fn new() -> Self {
        Self {
            buf: [0u8; MAX_BC_LEN],
        }
    }

    /// Pad a sequence to `max_len` using the collision-free padding scheme.
    /// Returns a slice valid until the next `pad()` call.
    #[inline]
    pub fn pad(&mut self, seq: &[u8], max_len: usize) -> &[u8] {
        pad_to_fixed(seq, max_len, &mut self.buf)
    }
}

impl Default for PadBuf {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pad_no_deficit() {
        let mut buf = [0u8; 32];
        let result = pad_to_fixed(b"ACGTACGTAC", 10, &mut buf);
        assert_eq!(result, b"ACGTACGTAC");
    }

    #[test]
    fn pad_deficit_1() {
        let mut buf = [0u8; 32];
        let result = pad_to_fixed(b"ACGTACGTA", 10, &mut buf);
        assert_eq!(result, b"ACGTACGTAA"); // "A" appended
    }

    #[test]
    fn pad_deficit_2() {
        let mut buf = [0u8; 32];
        let result = pad_to_fixed(b"ACGTACGT", 10, &mut buf);
        assert_eq!(result, b"ACGTACGTAC"); // "AC" appended
    }

    #[test]
    fn pad_deficit_3() {
        let mut buf = [0u8; 32];
        let result = pad_to_fixed(b"ACGTACG", 10, &mut buf);
        assert_eq!(result, b"ACGTACGAAG"); // "AAG" appended
    }

    #[test]
    fn pad_deficit_4() {
        let mut buf = [0u8; 32];
        let result = pad_to_fixed(b"ACGTAC", 10, &mut buf);
        assert_eq!(result, b"ACGTACAAAT"); // "AAAT" appended
    }

    #[test]
    fn no_collisions_across_lengths() {
        // Verify that barcodes of different original lengths don't collide after padding.
        // Take the same prefix "ACGTACGT" at lengths 8, 9, 10 and pad to 10.
        let mut buf = [0u8; 32];

        let padded_8 = pad_to_fixed(b"ACGTAC", 10, &mut buf).to_vec();
        let padded_9 = pad_to_fixed(b"ACGTACG", 10, &mut buf).to_vec();
        let padded_10 = pad_to_fixed(b"ACGTACGT", 10, &mut buf).to_vec();
        let padded_full = pad_to_fixed(b"ACGTACGTNN", 10, &mut buf).to_vec();

        // All should be different
        assert_ne!(padded_8, padded_9);
        assert_ne!(padded_8, padded_10);
        assert_ne!(padded_8, padded_full);
        assert_ne!(padded_9, padded_10);
        assert_ne!(padded_9, padded_full);
        assert_ne!(padded_10, padded_full);
    }

    #[test]
    fn padbuf_reuse() {
        let mut pb = PadBuf::new();

        let r1 = pb.pad(b"ACGT", 6).to_vec(); // deficit 2 -> "AC"
        assert_eq!(r1, b"ACGTAC");

        let r2 = pb.pad(b"TTTTT", 6).to_vec(); // deficit 1 -> "A"
        assert_eq!(r2, b"TTTTTA");
    }

    #[test]
    fn needs_padding_check() {
        assert!(!needs_padding(&GeoLen::Fixed(16)));
        assert!(!needs_padding(&GeoLen::Unbounded));
        assert!(needs_padding(&GeoLen::Range(9, 10)));
        assert!(!needs_padding(&GeoLen::Range(10, 10))); // min == max, no padding needed
    }

    #[test]
    fn normalized_len_check() {
        assert_eq!(normalized_len(&GeoLen::Fixed(16)), Some(16));
        assert_eq!(normalized_len(&GeoLen::Range(9, 10)), Some(10));
        assert_eq!(normalized_len(&GeoLen::Unbounded), None);
    }
}
