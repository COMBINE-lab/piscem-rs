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

/// The executor complexity tier required to interpret a geometry.
///
/// This separates geometries that can be handled by the current mostly
/// left-to-right executor from those that require a more general boundary
/// solver.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GeometryComplexity {
    /// Every field has a fixed width and can be extracted by static offsets.
    FixedOffsets,
    /// Exactly one variable-width region per read, and each such region is
    /// inferable from a fixed right boundary or the end of the read.
    ///
    /// This is the tier handled by the current executor.
    InferableVariable,
    /// A more general boundary-solving geometry, such as an interior `r:` or
    /// a variable region that must be inferred from both left and right
    /// boundaries after anchor resolution.
    ///
    /// Example: `1{r:f[ACAGT]b[9-11]}`.
    BoundarySolved,
}

/// A resolved boundary that can constrain a variable-width region.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BoundaryConstraint {
    /// The start of the read.
    ReadStart,
    /// The end of the read.
    ReadEnd,
    /// A fixed anchor sequence located within the read.
    Anchor(AnchorConstraint),
}

/// A fixed-sequence anchor used as a boundary.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AnchorConstraint {
    pub sequence: Vec<u8>,
    pub tolerance: Option<MatchTolerance>,
}

/// A bounded variable-width region inferred from explicit boundaries.
///
/// This is the natural unit for the current executor tier: once the left and
/// right boundaries are resolved, the width of the region is uniquely
/// determined and its inner tags can be sliced directly.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InferableRegion {
    pub left_boundary: BoundaryConstraint,
    pub right_boundary: BoundaryConstraint,
    pub parts: Vec<GeoPart>,
}

/// Sketch of a more general geometry representation for a boundary-solving
/// executor.
///
/// The intended model is that a read is decomposed into segments separated by
/// constraints, and extraction becomes a process of resolving boundary
/// positions and then assigning the spans between them.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BoundarySolvedReadGeom {
    pub segments: Vec<BoundarySolvedSegment>,
}

/// One segment in a boundary-solving read geometry.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BoundarySolvedSegment {
    /// A field or group of fields whose span is determined after solving its
    /// surrounding boundaries.
    Region(InferableRegion),
    /// A field that semantically consumes the maximal interval consistent with
    /// the remaining constraints, such as an interior `r:`.
    OpenEnded {
        tag: GeoTagType,
        left_boundary: BoundaryConstraint,
        right_boundary: BoundaryConstraint,
    },
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
