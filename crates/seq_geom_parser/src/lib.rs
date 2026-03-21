//! # seq_geom_parser
//!
//! A parser and extractor for sequencing read geometry descriptions.
//!
//! Geometry strings describe the layout of technical and biological sequences
//! within sequencing reads. For example, `1{b[16]u[12]x:}2{r:}` means:
//! - Read 1: 16bp cell barcode, 12bp UMI, discard rest
//! - Read 2: biological read (full length)
//!
//! ## Supported tags
//!
//! | Tag | Meaning | Example |
//! |-----|---------|---------|
//! | `b[N]` | Cell barcode | `b[16]` |
//! | `bN[L]` | Numbered barcode at level N | `b0[8]` |
//! | `s[N]` | Sample barcode (sugar for b0) | `s[8]` |
//! | `u[N]` | UMI | `u[12]` |
//! | `r[N]` / `r:` | Biological read (fixed/unbounded) | `r[50]`, `r:` |
//! | `f[SEQ]` | Fixed anchor sequence | `f[TTGCTAGGACCG]` |
//! | `x[N]` / `x:` | Discard (fixed/unbounded) | `x[18]`, `x:` |
//! | `x[N-M]` | Variable-length discard | `x[0-3]` |
//!
//! ## Distance functions
//!
//! Fixed anchors can be wrapped in distance functions for approximate matching:
//! - `hamming(f[SEQ], N)` — match within Hamming distance N
//!
//! ## Variable-length normalization
//!
//! Variable-length tags are supported when a downstream fixed anchor makes
//! their boundaries inferable, such as `b[9-10]u[12]f[SEQ]` or
//! `x[0-3]f[SEQ]s[10]`. Normalization helpers are exposed in [`normalize`] so
//! callers can pad extracted variable-length barcode/UMI sequences to their
//! declared maximum width when needed.
//!
//! ## Complexity Tiers
//!
//! The public API distinguishes three extraction tiers:
//! - [`GeometryComplexity::FixedOffsets`]: every extracted field has a static
//!   offset. Example: `1{b[16]u[12]x:}2{r:}`.
//! - [`GeometryComplexity::InferableVariable`]: one variable-width region per
//!   read, inferred from a fixed right boundary. Example:
//!   `1{b[9-10]f[ACGT]u[12]}2{r:}`.
//! - [`GeometryComplexity::BoundaryResolved`]: the read must first be split by
//!   resolved boundaries such as anchors and read ends. Example:
//!   `1{r:f[ACAGT]b[9-11]}2{u[12]x:}`.
//!
//! ## Boundary Resolution
//!
//! For [`GeometryComplexity::BoundaryResolved`] geometries, extraction proceeds
//! in two phases:
//! 1. Resolve anchor positions in read order.
//! 2. Assign the spans between those resolved boundaries to fields.
//!
//! If multiple anchor placements satisfy the geometry, the solver chooses the
//! monotone placement chain with the minimum total distance score. Ties are
//! broken by choosing the lexicographically leftmost anchor positions.
//!
//! ## Public Model vs Compiled Executor
//!
//! The boundary-oriented types exposed from [`types`] describe the public
//! semantic model of a geometry: boundaries, anchors, and segments between
//! resolved boundaries.
//!
//! They are not the same as the extractor's internal compiled representation.
//! [`CompiledGeom`] compiles parsed geometries into private extraction plans in
//! [`extract`] that are optimized for the hot path. This split is intentional:
//! the public types document the model and complexity hierarchy, while the
//! executor keeps a separate IR that can evolve for performance without
//! changing the public API.
//!
//! ## Examples By Tier
//!
//! ```rust
//! use seq_geom_parser::{geometry_complexity, parse_geometry, GeometryComplexity};
//!
//! let simple = parse_geometry("1{b[16]u[12]x:}2{r:}").unwrap();
//! assert_eq!(geometry_complexity(&simple), GeometryComplexity::FixedOffsets);
//!
//! let inferable = parse_geometry("1{b[9-10]f[ACGT]u[12]}2{r:}").unwrap();
//! assert_eq!(
//!     geometry_complexity(&inferable),
//!     GeometryComplexity::InferableVariable
//! );
//!
//! let boundary = parse_geometry("1{r:f[ACAGT]b[9-11]}2{u[12]x:}").unwrap();
//! assert_eq!(
//!     geometry_complexity(&boundary),
//!     GeometryComplexity::BoundaryResolved
//! );
//! ```

pub mod extract;
pub mod normalize;
pub mod parse;
pub mod types;

// Re-export key types at crate root
pub use extract::{
    BoundaryResolvedExtractor, CompiledGeom, ExtractedSeqs, GeomMeta, InferableExtractor,
    SimpleExtractor,
};
pub use parse::{format_errors, geometry_complexity, parse_geometry, validate_geometry};
pub use types::*;
