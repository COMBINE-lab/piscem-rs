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
//! The public types distinguish between two executor tiers:
//! - [`GeometryComplexity::FixedOffsets`]: fully static offsets
//! - [`GeometryComplexity::InferableVariable`]: one inferable variable region
//!   per read, resolved against an anchor or read boundary
//!
//! The crate also exposes sketch types for a future
//! [`GeometryComplexity::BoundarySolved`] executor, which would be needed for
//! broader grammars such as `1{r:f[ACAGT]b[9-11]}`.

pub mod extract;
pub mod normalize;
pub mod parse;
pub mod types;

// Re-export key types at crate root
pub use extract::{CompiledGeom, ExtractedSeqs, GeneralExtractor, GeomMeta, SimpleExtractor};
pub use parse::{format_errors, geometry_complexity, parse_geometry, validate_geometry};
pub use types::*;
