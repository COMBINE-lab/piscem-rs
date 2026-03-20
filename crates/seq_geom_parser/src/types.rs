//! Types for describing sequencing read geometry.

use smallvec::SmallVec;

/// Maximum number of barcode levels that can be stored inline (no heap allocation).
pub const MAX_INLINE_BARCODES: usize = 4;

/// The type/role of a geometry tag.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GeoTagType {
    /// Cell barcode (single-barcode protocols): `b[N]`
    Barcode,
    /// Numbered barcode at a specific level: `b0[N]`, `b1[N]`, etc.
    NumberedBarcode(u8),
    /// Sample/probe barcode (syntactic sugar for b0): `s[N]`
    SampleBarcode,
    /// Unique molecular identifier: `u[N]`
    Umi,
    /// Biological read sequence: `r[N]` or `r:`
    Read,
    /// Fixed/anchor DNA sequence: `f[ACGT...]`
    Fixed,
    /// Discard (skip) bases: `x[N]` or `x:`
    Discard,
}

/// The kind of distance metric for approximate matching.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DistanceKind {
    /// Hamming (substitution-only) distance.
    Hamming,
    /// Levenshtein (edit) distance — reserved for future use.
    Levenshtein,
}

/// Tolerance specification for approximate matching of fixed sequences.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MatchTolerance {
    pub kind: DistanceKind,
    pub max_dist: u8,
}

/// The length specification for a geometry tag.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GeoLen {
    /// Fixed length: `[N]`
    Fixed(u32),
    /// Variable-length range: `[N-M]` where N <= M
    Range(u32, u32),
    /// Unbounded (rest of read): `:`
    Unbounded,
}

/// A single piece of a geometry description.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GeoPart {
    pub tag: GeoTagType,
    pub len: GeoLen,
    /// For `Fixed` tags: the expected DNA sequence (uppercase ACGT).
    pub sequence: Option<Vec<u8>>,
    /// For `Fixed` tags wrapped in `hamming(...)` etc.: the matching tolerance.
    pub tolerance: Option<MatchTolerance>,
}

/// Describes the geometry of one read (R1 or R2).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReadGeom {
    pub parts: Vec<GeoPart>,
}

/// A fully parsed geometry description for a paired-end library.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FragmentGeom {
    pub read1: ReadGeom,
    pub read2: ReadGeom,
}

/// Information about barcode levels discovered in the geometry.
#[derive(Debug, Clone)]
pub struct BarcodeInfo {
    /// Number of barcode levels (1 for standard, 2+ for multi-barcode).
    pub num_levels: usize,
    /// The role of each barcode level.
    pub roles: SmallVec<[BarcodeRole; MAX_INLINE_BARCODES]>,
}

/// The semantic role of a barcode level.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BarcodeRole {
    /// Sample/library barcode (from `s[]` or explicit `b0[]` with sample semantics).
    Sample,
    /// Cell barcode.
    Cell,
    /// Generic numbered barcode with no specific semantic role.
    Generic(u8),
}

/// Padding table for variable-length barcode normalization.
/// Sequences of length `max - k` (for k = 0..4) are padded with these suffixes.
/// The suffixes are chosen so that no padded barcode of length L collides with
/// a padded barcode of length L' when L != L'.
pub const VAR_LEN_PADDING: &[&[u8]] = &[
    b"",     // captured_len == max: no padding
    b"A",    // captured_len == max - 1
    b"AC",   // captured_len == max - 2
    b"AAG",  // captured_len == max - 3
    b"AAAT", // captured_len == max - 4
];

/// Maximum allowed range width (max - min) for variable-length tags.
pub const MAX_RANGE_WIDTH: u32 = 4;
