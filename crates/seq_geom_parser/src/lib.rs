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
//! When barcodes or UMIs have variable length (`b[9-10]`, `u[10-12]`), the
//! extracted sequence is padded to the maximum length using a collision-free
//! padding scheme, ensuring downstream fixed-width tools work correctly.

pub mod types;
pub mod parse;
pub mod extract;
pub mod normalize;

// Re-export key types at crate root
pub use types::*;
pub use parse::{parse_geometry, validate_geometry, format_errors};
pub use extract::{CompiledGeom, SimpleExtractor, GeneralExtractor, GeomMeta, ExtractedSeqs};
