//! Sequence extraction from reads using a compiled geometry.
//!
//! The extraction is split into two phases:
//! 1. **Compile**: Convert a parsed [`FragmentGeom`] into a [`CompiledGeom`] with
//!    precomputed offsets and extraction plans. This happens once at startup.
//! 2. **Extract**: Apply the compiled geometry to each read pair. This is the hot
//!    path and must be as fast as possible — zero allocation, minimal branching.
//!
//! For fixed-length geometries (the common case), extraction is a simple slice
//! operation: `&read[offset..offset+len]`. For inferable variable-length
//! geometries, a search phase locates the right boundary (an anchor or the read
//! end) and resolves the variable-width span before slicing.
#![allow(clippy::too_many_arguments, clippy::type_complexity)]

use smallvec::{smallvec, SmallVec};

use crate::parse::geometry_complexity;
use crate::types::*;

/// Normalization requirements for extracted technical sequences.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NormalizationMeta {
    /// Whether each barcode level may require padding.
    pub bc_needs_padding: SmallVec<[bool; MAX_INLINE_BARCODES]>,
    /// Whether the UMI may require padding.
    pub umi_needs_padding: bool,
    /// Whether any barcode level or the UMI may require padding.
    pub any_normalization: bool,
}

/// A precomputed extraction plan for one read.
#[derive(Debug, Clone)]
struct ReadPlan {
    /// Extraction steps, in order. Each step extracts or skips a region.
    steps: Vec<ExtractionStep>,
    /// Whether this read contains any variable-length steps (requires search).
    has_variable: bool,
}

#[derive(Debug, Clone)]
struct BoundaryResolvedReadPlan {
    anchors: Vec<AnchorPlan>,
    segments: Vec<ResolvedSegmentPlan>,
}

#[derive(Debug, Clone)]
struct AnchorPlan {
    sequence: Vec<u8>,
    max_dist: u8,
    dist_kind: DistanceKind,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SegmentBoundaryRef {
    ReadStart,
    ReadEnd,
    AnchorStart(usize),
    AnchorEnd(usize),
}

#[derive(Debug, Clone)]
struct ResolvedSegmentPlan {
    left: SegmentBoundaryRef,
    right: SegmentBoundaryRef,
    parts: Vec<ResolvedSegmentPart>,
}

#[derive(Debug, Clone)]
enum ResolvedSegmentPart {
    Fixed {
        len: usize,
        target: Option<ExtractTarget>,
    },
    Variable {
        min_len: usize,
        max_len: usize,
        target: Option<ExtractTarget>,
    },
    Unbounded {
        target: Option<ExtractTarget>,
    },
}

/// A single extraction step within a read.
#[derive(Debug, Clone)]
enum ExtractionStep {
    /// Extract a fixed-offset slice: `&read[offset..offset+len]`
    FixedSlice {
        offset: usize,
        len: usize,
        target: ExtractTarget,
    },
    /// Skip a fixed number of bases (no extraction).
    FixedSkip,
    /// Extract from offset to end of read.
    Unbounded {
        offset: usize,
        target: ExtractTarget,
    },
    /// Skip to end of read (unbounded discard).
    UnboundedSkip,
    /// Search for an anchor sequence within a window, then continue
    /// extraction relative to the anchor position.
    AnchorSearch {
        /// Absolute offset where the bounded region starts.
        region_start_offset: usize,
        /// Minimum offset to start searching.
        min_offset: usize,
        /// Maximum offset to start searching (inclusive).
        max_offset: usize,
        /// The anchor sequence to search for.
        anchor: Vec<u8>,
        /// Maximum allowed distance for approximate matching (0 = exact).
        max_dist: u8,
        /// Distance metric kind.
        dist_kind: DistanceKind,
        /// Steps within the bounded region immediately preceding the anchor.
        pre_anchor_parts: Vec<AnchoredPart>,
        /// Steps to execute after the anchor is found. Offsets are relative
        /// to the position immediately after the anchor.
        post_anchor_steps: Vec<ExtractionStep>,
    },
}

#[derive(Debug, Clone)]
enum AnchoredPart {
    Fixed {
        len: usize,
        target: Option<ExtractTarget>,
    },
    Variable {
        min_len: usize,
        max_len: usize,
        target: Option<ExtractTarget>,
    },
}

/// What to extract a slice into.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ExtractTarget {
    /// Barcode at level N.
    Barcode(u8),
    /// UMI.
    Umi,
    /// Biological read.
    Read,
}

/// Result of extracting sequences from a read pair.
#[derive(Debug)]
pub struct ExtractedSeqs<'a> {
    /// Barcodes by level. For single-barcode protocols, only index 0 is used.
    pub barcodes: SmallVec<[Option<&'a [u8]>; MAX_INLINE_BARCODES]>,
    /// UMI sequence.
    pub umi: Option<&'a [u8]>,
    /// Biological read(s). For SE: one entry. For PE: two entries.
    pub reads: SmallVec<[&'a [u8]; 2]>,
}

// NOTE ON EXTRACTION API DESIGN
//
// We benchmarked several approaches for the per-read extraction hot path:
//
//   Hardcoded (tuple return):    1.37 ns  — absolute minimum, reference only
//   simple.extract() by-value:   3.18 ns  — FASTEST generic approach
//   enum dispatch:               3.57 ns  — convenience, negligible overhead
//   extract_into(&mut buf):      4.60 ns  — reusable buffer, slower (pointer indirection)
//   Stateful extractor + getters: 6.70 ns — raw pointer storage, worst of all
//
// The by-value return wins because for 1-4 barcode levels, SmallVec stores
// data inline (on the stack), and the compiler optimizes the construction
// into register operations. Any attempt to reuse storage via mutable
// references or raw pointers introduces memory indirection that the
// optimizer can't eliminate, making it slower than just constructing a
// fresh SmallVec each time.
//
// The recommended usage pattern:
//
//   match CompiledGeom::from_fragment_geom(&geom)? {
//       CompiledGeom::Simple(ext) => {
//           for (r1, r2) in reads {
//               let seqs = ext.extract(r1, r2); // 3.18 ns, no dispatch
//           }
//       }
//       CompiledGeom::Inferable(ext) => {
//           for (r1, r2) in reads {
//               let seqs = ext.extract(r1, r2); // general path
//           }
//       }
//   }
//

/// Specialized fast extraction for the common case:
/// R1 has BC at [0..bc_len], UMI at [bc_len..bc_len+umi_len], R2 is the full bio read.
/// This avoids all step iteration, SmallVec construction, and enum matching.
#[derive(Debug, Clone, Copy)]
struct SimpleGeom {
    bc_offset: usize,
    bc_len: usize,
    umi_offset: usize,
    umi_len: usize,
    /// true = bio read on R2, false = bio read on R1 (after tech prefix)
    read_on_r2: bool,
    /// If bio read is on R1, the start offset
    read_offset: usize,
}

/// A compiled geometry ready for fast extraction.
///
/// This is an enum so that callers can match once at startup and then call
/// the variant-specific `extract` method in a tight loop with no per-read
/// branching. The branch predictor would handle the enum match on `extract()`
/// perfectly (variant never changes), but callers who want guaranteed zero
/// overhead can match once and hold the inner type directly.
///
/// Variant selection:
/// - [`CompiledGeom::Simple`] for fixed-offset layouts such as
///   `1{b[16]u[12]x:}2{r:}`
/// - [`CompiledGeom::Inferable`] for locally bounded variable-width layouts
///   such as `1{b[9-10]f[ACGT]u[12]}2{r:}`
/// - [`CompiledGeom::BoundaryResolved`] for layouts that require anchor
///   placement before segments can be assigned, such as
///   `1{r:f[ACAGT]b[9-11]}2{u[12]x:}`
#[derive(Debug, Clone)]
pub enum CompiledGeom {
    /// Fast path for simple single-barcode geometries (e.g. Chromium v2/v3).
    /// BC+UMI on one read at fixed offsets, bio read on the other.
    Simple(SimpleExtractor),
    /// Inferable path for anchored or otherwise locally bounded variable-width geometries.
    Inferable(InferableExtractor),
    /// Boundary-resolved path for geometries that require anchor placement
    /// before open-ended or variable-width segments can be assigned.
    BoundaryResolved(BoundaryResolvedExtractor),
}

/// Metadata common to all extractor variants.
#[derive(Debug, Clone)]
pub struct GeomMeta {
    /// Number of barcode levels.
    pub num_bc_levels: usize,
    /// Barcode lengths (fixed, after normalization). One per level.
    pub bc_lens: SmallVec<[usize; MAX_INLINE_BARCODES]>,
    /// UMI length (fixed, after normalization).
    pub umi_len: usize,
    /// Whether the biological read is paired.
    pub is_paired_read: bool,
    /// Whether any barcode/UMI needs padding before downstream fixed-width encoding.
    pub normalization: NormalizationMeta,
}

/// Fast extractor for simple single-barcode geometries.
#[derive(Debug, Clone)]
pub struct SimpleExtractor {
    pub meta: GeomMeta,
    sg: SimpleGeom,
}

/// Inferable extractor for anchored variable-width geometries.
///
/// This tier assumes a read can still be executed mostly left-to-right once a
/// single right boundary has been located.
#[derive(Debug, Clone)]
pub struct InferableExtractor {
    pub meta: GeomMeta,
    r1_plan: ReadPlan,
    r2_plan: ReadPlan,
}

/// Extractor for boundary-resolved geometries with interior open-ended spans.
///
/// The solver first resolves anchors in read order, choosing the monotone chain
/// with minimum total distance score. If several chains have the same score, it
/// chooses the lexicographically leftmost anchor positions. The spans between
/// those resolved boundaries are then assigned to the segment fields.
#[derive(Debug, Clone)]
pub struct BoundaryResolvedExtractor {
    pub meta: GeomMeta,
    r1_plan: BoundaryResolvedReadPlan,
    r2_plan: BoundaryResolvedReadPlan,
}

impl SimpleExtractor {
    /// Extract sequences from a read pair. Zero-cost: no branching, no step iteration.
    #[inline(always)]
    pub fn extract<'a>(&self, r1: &'a [u8], r2: &'a [u8]) -> ExtractedSeqs<'a> {
        let sg = &self.sg;
        let barcode = if r1.len() >= sg.bc_offset + sg.bc_len {
            Some(&r1[sg.bc_offset..sg.bc_offset + sg.bc_len])
        } else {
            None
        };
        let umi = if r1.len() >= sg.umi_offset + sg.umi_len {
            Some(&r1[sg.umi_offset..sg.umi_offset + sg.umi_len])
        } else {
            None
        };
        let mut reads = SmallVec::new();
        if sg.read_on_r2 {
            reads.push(r2);
        } else if sg.read_offset < r1.len() {
            reads.push(&r1[sg.read_offset..]);
        }
        ExtractedSeqs {
            barcodes: smallvec![barcode],
            umi,
            reads,
        }
    }
}

impl InferableExtractor {
    /// Extract sequences from a read pair using step-based plans.
    #[inline]
    pub fn extract<'a>(&self, r1: &'a [u8], r2: &'a [u8]) -> ExtractedSeqs<'a> {
        let mut result = ExtractedSeqs {
            barcodes: SmallVec::from_elem(None, self.meta.num_bc_levels),
            umi: None,
            reads: SmallVec::new(),
        };
        execute_plan(&self.r1_plan, r1, &mut result);
        execute_plan(&self.r2_plan, r2, &mut result);
        result
    }
}

impl BoundaryResolvedExtractor {
    #[inline]
    pub fn extract<'a>(&self, r1: &'a [u8], r2: &'a [u8]) -> ExtractedSeqs<'a> {
        let mut result = ExtractedSeqs {
            barcodes: SmallVec::from_elem(None, self.meta.num_bc_levels),
            umi: None,
            reads: SmallVec::new(),
        };
        execute_boundary_resolved_plan(&self.r1_plan, r1, &mut result);
        execute_boundary_resolved_plan(&self.r2_plan, r2, &mut result);
        result
    }
}

fn build_geom_meta(
    r1_bc_info: &[usize],
    r2_bc_info: &[usize],
    r1_bc_padding: &[bool],
    r2_bc_padding: &[bool],
    r1_umi_len: usize,
    r2_umi_len: usize,
    r1_umi_needs_padding: bool,
    r2_umi_needs_padding: bool,
    r1_has_read: bool,
    r2_has_read: bool,
) -> GeomMeta {
    let mut bc_lens = SmallVec::new();
    let mut bc_needs_padding = SmallVec::new();
    let num_bc_levels = r1_bc_info.len().max(r2_bc_info.len());
    for level in 0..num_bc_levels {
        let r1_len = r1_bc_info.get(level).copied().unwrap_or(0);
        let r2_len = r2_bc_info.get(level).copied().unwrap_or(0);
        bc_lens.push(r1_len.max(r2_len));

        let r1_pad = r1_bc_padding.get(level).copied().unwrap_or(false);
        let r2_pad = r2_bc_padding.get(level).copied().unwrap_or(false);
        bc_needs_padding.push(r1_pad || r2_pad);
    }

    let umi_needs_padding = r1_umi_needs_padding || r2_umi_needs_padding;
    let any_normalization = umi_needs_padding || bc_needs_padding.iter().copied().any(|v| v);

    GeomMeta {
        num_bc_levels,
        bc_lens,
        umi_len: r1_umi_len.max(r2_umi_len),
        is_paired_read: r1_has_read && r2_has_read,
        normalization: NormalizationMeta {
            bc_needs_padding,
            umi_needs_padding,
            any_normalization,
        },
    }
}

fn detect_simple_geom(r1_plan: &ReadPlan, r2_plan: &ReadPlan) -> Option<SimpleGeom> {
    let mut bc_info = None;
    let mut umi_info = None;
    let mut r1_read_offset = None;
    for step in &r1_plan.steps {
        match step {
            ExtractionStep::FixedSlice {
                offset,
                len,
                target: ExtractTarget::Barcode(0),
            } => bc_info = Some((*offset, *len)),
            ExtractionStep::FixedSlice {
                offset,
                len,
                target: ExtractTarget::Umi,
            } => umi_info = Some((*offset, *len)),
            ExtractionStep::Unbounded {
                offset,
                target: ExtractTarget::Read,
            } => r1_read_offset = Some(*offset),
            _ => {}
        }
    }

    let r2_has_unbounded_read = r2_plan.steps.iter().any(|s| {
        matches!(
            s,
            ExtractionStep::Unbounded {
                offset: 0,
                target: ExtractTarget::Read
            }
        )
    });

    match (bc_info, umi_info) {
        (Some((bc_off, bc_l)), Some((umi_off, umi_l))) if r2_has_unbounded_read => {
            Some(SimpleGeom {
                bc_offset: bc_off,
                bc_len: bc_l,
                umi_offset: umi_off,
                umi_len: umi_l,
                read_on_r2: true,
                read_offset: 0,
            })
        }
        (Some((bc_off, bc_l)), Some((umi_off, umi_l))) if r1_read_offset.is_some() => {
            Some(SimpleGeom {
                bc_offset: bc_off,
                bc_len: bc_l,
                umi_offset: umi_off,
                umi_len: umi_l,
                read_on_r2: false,
                read_offset: r1_read_offset.unwrap(),
            })
        }
        _ => None,
    }
}

impl CompiledGeom {
    /// Compile a parsed geometry into an extraction plan.
    pub fn from_fragment_geom(geom: &FragmentGeom) -> Result<Self, String> {
        let bc_level_fn = compute_barcode_levels(geom);
        match geometry_complexity(geom) {
            GeometryComplexity::FixedOffsets | GeometryComplexity::InferableVariable => {
                let (
                    r1_plan,
                    r1_bc_info,
                    r1_bc_padding,
                    r1_umi_len,
                    r1_umi_needs_padding,
                    r1_has_read,
                ) = compile_read_plan_with_levels(&geom.read1, &bc_level_fn)?;
                let (
                    r2_plan,
                    r2_bc_info,
                    r2_bc_padding,
                    r2_umi_len,
                    r2_umi_needs_padding,
                    r2_has_read,
                ) = compile_read_plan_with_levels(&geom.read2, &bc_level_fn)?;

                let meta = build_geom_meta(
                    &r1_bc_info,
                    &r2_bc_info,
                    &r1_bc_padding,
                    &r2_bc_padding,
                    r1_umi_len,
                    r2_umi_len,
                    r1_umi_needs_padding,
                    r2_umi_needs_padding,
                    r1_has_read,
                    r2_has_read,
                );

                let simple =
                    if meta.num_bc_levels == 1 && !r1_plan.has_variable && !r2_plan.has_variable {
                        detect_simple_geom(&r1_plan, &r2_plan)
                    } else {
                        None
                    };

                if let Some(sg) = simple {
                    Ok(CompiledGeom::Simple(SimpleExtractor { meta, sg }))
                } else {
                    Ok(CompiledGeom::Inferable(InferableExtractor {
                        meta,
                        r1_plan,
                        r2_plan,
                    }))
                }
            }
            GeometryComplexity::BoundaryResolved => {
                let (
                    r1_plan,
                    r1_bc_info,
                    r1_bc_padding,
                    r1_umi_len,
                    r1_umi_needs_padding,
                    r1_has_read,
                ) = compile_boundary_resolved_read_plan(&geom.read1, &bc_level_fn)?;
                let (
                    r2_plan,
                    r2_bc_info,
                    r2_bc_padding,
                    r2_umi_len,
                    r2_umi_needs_padding,
                    r2_has_read,
                ) = compile_boundary_resolved_read_plan(&geom.read2, &bc_level_fn)?;

                let meta = build_geom_meta(
                    &r1_bc_info,
                    &r2_bc_info,
                    &r1_bc_padding,
                    &r2_bc_padding,
                    r1_umi_len,
                    r2_umi_len,
                    r1_umi_needs_padding,
                    r2_umi_needs_padding,
                    r1_has_read,
                    r2_has_read,
                );

                Ok(CompiledGeom::BoundaryResolved(BoundaryResolvedExtractor {
                    meta,
                    r1_plan,
                    r2_plan,
                }))
            }
        }
    }

    /// Extract sequences from a read pair.
    ///
    /// This dispatches to the variant-specific extractor. For maximum
    /// performance in tight loops, match on the enum once at startup and
    /// call the variant's `extract` method directly — this eliminates even
    /// the (perfectly-predicted) enum branch.
    ///
    /// ```ignore
    /// match compiled {
    ///     CompiledGeom::Simple(ext) => {
    ///         for (r1, r2) in reads {
    ///             let seqs = ext.extract(r1, r2); // no dispatch
    ///         }
    ///     }
    ///     CompiledGeom::Inferable(ext) => {
    ///         for (r1, r2) in reads {
    ///             let seqs = ext.extract(r1, r2);
    ///         }
    ///     }
    ///     CompiledGeom::BoundaryResolved(ext) => {
    ///         for (r1, r2) in reads {
    ///             let seqs = ext.extract(r1, r2);
    ///         }
    ///     }
    /// }
    /// ```
    #[inline]
    pub fn extract<'a>(&self, r1: &'a [u8], r2: &'a [u8]) -> ExtractedSeqs<'a> {
        match self {
            CompiledGeom::Simple(ext) => ext.extract(r1, r2),
            CompiledGeom::Inferable(ext) => ext.extract(r1, r2),
            CompiledGeom::BoundaryResolved(ext) => ext.extract(r1, r2),
        }
    }

    /// Access geometry metadata (barcode lengths, UMI length, etc.).
    pub fn meta(&self) -> &GeomMeta {
        match self {
            CompiledGeom::Simple(ext) => &ext.meta,
            CompiledGeom::Inferable(ext) => &ext.meta,
            CompiledGeom::BoundaryResolved(ext) => &ext.meta,
        }
    }
}

/// Execute an extraction plan against a read, filling in the result.
#[inline]
fn execute_plan<'a>(plan: &ReadPlan, read: &'a [u8], result: &mut ExtractedSeqs<'a>) {
    if !plan.has_variable {
        // Fast path: all steps are fixed-offset. No searching needed.
        for step in &plan.steps {
            execute_fixed_step(step, read, result);
        }
    } else {
        // Slow path: has variable-length steps with anchor search.
        execute_steps_with_search(&plan.steps, read, result);
    }
}

fn execute_boundary_resolved_plan<'a>(
    plan: &BoundaryResolvedReadPlan,
    read: &'a [u8],
    result: &mut ExtractedSeqs<'a>,
) {
    let Some(anchor_positions) = resolve_anchor_positions(plan, read) else {
        return;
    };

    for segment in &plan.segments {
        execute_resolved_segment(segment, &anchor_positions, read, result);
    }
}

/// Execute a single fixed-offset step (no branching on variable-length).
#[inline(always)]
fn execute_fixed_step<'a>(step: &ExtractionStep, read: &'a [u8], result: &mut ExtractedSeqs<'a>) {
    match step {
        ExtractionStep::FixedSlice {
            offset,
            len,
            target,
        } => {
            let end = *offset + *len;
            if end <= read.len() {
                let slice = &read[*offset..end];
                assign_target(target, slice, result);
            }
        }
        ExtractionStep::Unbounded { offset, target } => {
            if *offset < read.len() {
                assign_target(target, &read[*offset..], result);
            }
        }
        ExtractionStep::FixedSkip | ExtractionStep::UnboundedSkip => {
            // Nothing to extract
        }
        ExtractionStep::AnchorSearch { .. } => {
            // Should not appear in fixed-only plans
        }
    }
}

/// Execute steps that may include anchor searches.
fn execute_steps_with_search<'a>(
    steps: &[ExtractionStep],
    read: &'a [u8],
    result: &mut ExtractedSeqs<'a>,
) {
    for step in steps {
        match step {
            ExtractionStep::AnchorSearch {
                region_start_offset,
                min_offset,
                max_offset,
                anchor,
                max_dist,
                dist_kind,
                pre_anchor_parts,
                post_anchor_steps,
            } => {
                // Search for the anchor within the window
                if let Some(anchor_pos) = find_anchor(
                    read,
                    *min_offset,
                    *max_offset,
                    anchor,
                    *max_dist,
                    *dist_kind,
                ) {
                    execute_anchored_parts(
                        *region_start_offset,
                        anchor_pos,
                        pre_anchor_parts,
                        read,
                        result,
                    );
                    let after_anchor = anchor_pos + anchor.len();
                    // Execute post-anchor steps with offsets relative to after_anchor
                    for post_step in post_anchor_steps {
                        execute_fixed_step_with_base(post_step, read, after_anchor, result);
                    }
                }
                // If anchor not found, skip all post-anchor extractions (read is unmatchable)
            }
            other => execute_fixed_step(other, read, result),
        }
    }
}

/// Execute a fixed step with a base offset added to all positions.
#[inline(always)]
fn execute_fixed_step_with_base<'a>(
    step: &ExtractionStep,
    read: &'a [u8],
    base: usize,
    result: &mut ExtractedSeqs<'a>,
) {
    match step {
        ExtractionStep::FixedSlice {
            offset,
            len,
            target,
        } => {
            let abs_offset = base + *offset;
            let end = abs_offset + *len;
            if end <= read.len() {
                assign_target(target, &read[abs_offset..end], result);
            }
        }
        ExtractionStep::Unbounded { offset, target } => {
            let abs_offset = base + *offset;
            if abs_offset < read.len() {
                assign_target(target, &read[abs_offset..], result);
            }
        }
        ExtractionStep::FixedSkip | ExtractionStep::UnboundedSkip => {}
        ExtractionStep::AnchorSearch { .. } => {
            // Nested anchor searches not supported
        }
    }
}

fn execute_anchored_parts<'a>(
    region_start_offset: usize,
    anchor_pos: usize,
    parts: &[AnchoredPart],
    read: &'a [u8],
    result: &mut ExtractedSeqs<'a>,
) {
    let total_span = match anchor_pos.checked_sub(region_start_offset) {
        Some(span) => span,
        None => return,
    };

    let fixed_total = parts
        .iter()
        .map(|part| match part {
            AnchoredPart::Fixed { len, .. } => *len,
            AnchoredPart::Variable { .. } => 0,
        })
        .sum::<usize>();

    let variable_part = parts.iter().find_map(|part| match part {
        AnchoredPart::Variable {
            min_len, max_len, ..
        } => Some((*min_len, *max_len)),
        AnchoredPart::Fixed { .. } => None,
    });

    let variable_len = if let Some((min_len, max_len)) = variable_part {
        let len = match total_span.checked_sub(fixed_total) {
            Some(len) => len,
            None => return,
        };
        if len < min_len || len > max_len {
            return;
        }
        len
    } else if total_span == fixed_total {
        0
    } else {
        return;
    };

    let mut offset = region_start_offset;
    for part in parts {
        match part {
            AnchoredPart::Fixed { len, target } => {
                let end = offset + *len;
                if end > read.len() {
                    return;
                }
                if let Some(target) = target {
                    assign_target(target, &read[offset..end], result);
                }
                offset = end;
            }
            AnchoredPart::Variable { target, .. } => {
                let end = offset + variable_len;
                if end > read.len() {
                    return;
                }
                if let Some(target) = target {
                    assign_target(target, &read[offset..end], result);
                }
                offset = end;
            }
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct AnchorCandidate {
    start: usize,
    end: usize,
    score: u8,
}

type AnchorChainState = Option<(u32, Vec<AnchorCandidate>)>;

fn execute_resolved_segment<'a>(
    segment: &ResolvedSegmentPlan,
    anchors: &[AnchorCandidate],
    read: &'a [u8],
    result: &mut ExtractedSeqs<'a>,
) {
    let Some(left) = boundary_position(segment.left, anchors, false, read.len()) else {
        return;
    };
    let Some(right) = boundary_position(segment.right, anchors, true, read.len()) else {
        return;
    };
    if right < left {
        return;
    }

    let span = right - left;
    let fixed_total = segment
        .parts
        .iter()
        .map(|part| match part {
            ResolvedSegmentPart::Fixed { len, .. } => *len,
            ResolvedSegmentPart::Variable { .. } | ResolvedSegmentPart::Unbounded { .. } => 0,
        })
        .sum::<usize>();

    let flexible = segment.parts.iter().find_map(|part| match part {
        ResolvedSegmentPart::Fixed { .. } => None,
        ResolvedSegmentPart::Variable {
            min_len, max_len, ..
        } => Some((*min_len, Some(*max_len))),
        ResolvedSegmentPart::Unbounded { .. } => Some((0usize, None)),
    });

    let flex_len = if let Some((min_len, max_len)) = flexible {
        let len = match span.checked_sub(fixed_total) {
            Some(len) => len,
            None => return,
        };
        if len < min_len {
            return;
        }
        if let Some(max_len) = max_len {
            if len > max_len {
                return;
            }
        }
        len
    } else if span == fixed_total {
        0
    } else {
        return;
    };

    let mut offset = left;
    for part in &segment.parts {
        match part {
            ResolvedSegmentPart::Fixed { len, target } => {
                let end = offset + *len;
                if end > read.len() {
                    return;
                }
                if let Some(target) = target {
                    assign_target(target, &read[offset..end], result);
                }
                offset = end;
            }
            ResolvedSegmentPart::Variable { target, .. }
            | ResolvedSegmentPart::Unbounded { target } => {
                let end = offset + flex_len;
                if end > read.len() {
                    return;
                }
                if let Some(target) = target {
                    assign_target(target, &read[offset..end], result);
                }
                offset = end;
            }
        }
    }
}

fn boundary_position(
    boundary: SegmentBoundaryRef,
    anchors: &[AnchorCandidate],
    _is_right: bool,
    read_len: usize,
) -> Option<usize> {
    match boundary {
        SegmentBoundaryRef::ReadStart => Some(0),
        SegmentBoundaryRef::ReadEnd => Some(read_len),
        SegmentBoundaryRef::AnchorStart(idx) => anchors.get(idx).map(|a| a.start),
        SegmentBoundaryRef::AnchorEnd(idx) => anchors.get(idx).map(|a| a.end),
    }
}

fn resolve_anchor_positions(
    plan: &BoundaryResolvedReadPlan,
    read: &[u8],
) -> Option<Vec<AnchorCandidate>> {
    if plan.anchors.is_empty() {
        return Some(Vec::new());
    }

    let candidates = plan
        .anchors
        .iter()
        .map(|anchor| find_anchor_candidates(read, anchor))
        .collect::<Vec<_>>();
    if candidates.iter().any(|c| c.is_empty()) {
        return None;
    }

    let seg_bounds = plan
        .segments
        .iter()
        .map(segment_span_bounds)
        .collect::<Vec<_>>();

    let mut states: Vec<Vec<AnchorChainState>> = vec![Vec::new(); candidates.len()];
    for (idx, cand) in candidates[0].iter().enumerate() {
        let prefix_span = cand.start;
        if span_satisfies(prefix_span, seg_bounds[0]) {
            states[0].resize(candidates[0].len(), None);
            states[0][idx] = Some((cand.score as u32, vec![*cand]));
        }
    }

    for anchor_idx in 1..candidates.len() {
        states[anchor_idx].resize(candidates[anchor_idx].len(), None);
        for (cand_idx, cand) in candidates[anchor_idx].iter().enumerate() {
            let mut best: Option<(u32, Vec<AnchorCandidate>)> = None;
            for (prev_idx, prev) in candidates[anchor_idx - 1].iter().enumerate() {
                let Some((prev_score, prev_chain)) = &states[anchor_idx - 1][prev_idx] else {
                    continue;
                };
                if prev.end > cand.start {
                    continue;
                }
                let span = cand.start - prev.end;
                if !span_satisfies(span, seg_bounds[anchor_idx]) {
                    continue;
                }

                let mut chain = prev_chain.clone();
                chain.push(*cand);
                let score = *prev_score + cand.score as u32;
                if best.as_ref().is_none_or(|(best_score, best_chain)| {
                    score < *best_score
                        || (score == *best_score && chain_lex_less(&chain, best_chain))
                }) {
                    best = Some((score, chain));
                }
            }
            states[anchor_idx][cand_idx] = best;
        }
    }

    let tail_idx = seg_bounds.len() - 1;
    let mut best_final: Option<(u32, Vec<AnchorCandidate>)> = None;
    let last_anchor_idx = candidates.len() - 1;
    for (cand_idx, cand) in candidates[last_anchor_idx].iter().enumerate() {
        let Some((score, chain)) = &states[last_anchor_idx][cand_idx] else {
            continue;
        };
        let tail_span = read.len().saturating_sub(cand.end);
        if !span_satisfies(tail_span, seg_bounds[tail_idx]) {
            continue;
        }
        if best_final.as_ref().is_none_or(|(best_score, best_chain)| {
            *score < *best_score || (*score == *best_score && chain_lex_less(chain, best_chain))
        }) {
            best_final = Some((*score, chain.clone()));
        }
    }

    best_final.map(|(_, chain)| chain)
}

fn chain_lex_less(a: &[AnchorCandidate], b: &[AnchorCandidate]) -> bool {
    a.iter()
        .zip(b.iter())
        .find_map(|(a, b)| {
            if a.start == b.start {
                None
            } else {
                Some(a.start < b.start)
            }
        })
        .unwrap_or(a.len() < b.len())
}

fn segment_span_bounds(segment: &ResolvedSegmentPlan) -> (usize, Option<usize>) {
    let fixed_total = segment
        .parts
        .iter()
        .map(|part| match part {
            ResolvedSegmentPart::Fixed { len, .. } => *len,
            ResolvedSegmentPart::Variable { .. } | ResolvedSegmentPart::Unbounded { .. } => 0,
        })
        .sum::<usize>();
    let flexible = segment.parts.iter().find_map(|part| match part {
        ResolvedSegmentPart::Fixed { .. } => None,
        ResolvedSegmentPart::Variable {
            min_len, max_len, ..
        } => Some((*min_len, Some(*max_len))),
        ResolvedSegmentPart::Unbounded { .. } => Some((0usize, None)),
    });

    match flexible {
        Some((min, Some(max))) => (fixed_total + min, Some(fixed_total + max)),
        Some((min, None)) => (fixed_total + min, None),
        None => (fixed_total, Some(fixed_total)),
    }
}

fn span_satisfies(span: usize, bounds: (usize, Option<usize>)) -> bool {
    span >= bounds.0 && bounds.1.is_none_or(|max| span <= max)
}

#[inline(always)]
fn assign_target<'a>(target: &ExtractTarget, slice: &'a [u8], result: &mut ExtractedSeqs<'a>) {
    match target {
        ExtractTarget::Barcode(level) => {
            if (*level as usize) < result.barcodes.len() {
                result.barcodes[*level as usize] = Some(slice);
            }
        }
        ExtractTarget::Umi => {
            result.umi = Some(slice);
        }
        ExtractTarget::Read => {
            result.reads.push(slice);
        }
    }
}

/// Search for an anchor sequence within a window of the read.
///
/// Returns the starting position of the best match, or None if no match
/// is within the distance tolerance.
#[inline]
fn find_anchor(
    read: &[u8],
    min_offset: usize,
    max_offset: usize,
    anchor: &[u8],
    max_dist: u8,
    dist_kind: DistanceKind,
) -> Option<usize> {
    let anchor_len = anchor.len();

    if max_dist == 0 {
        // Exact match: simple byte comparison
        for pos in min_offset..=max_offset {
            if pos + anchor_len <= read.len() && read[pos..pos + anchor_len] == *anchor {
                return Some(pos);
            }
        }
        None
    } else {
        // Approximate match: find the position with minimum distance
        let mut best_pos = None;
        let mut best_dist = max_dist + 1;

        for pos in min_offset..=max_offset {
            if pos + anchor_len > read.len() {
                continue;
            }
            let d = match dist_kind {
                DistanceKind::Hamming => {
                    hamming_distance_at_most(&read[pos..pos + anchor_len], anchor, max_dist)
                }
                DistanceKind::Levenshtein => {
                    // TODO: implement when needed
                    u8::MAX
                }
            };
            if d < best_dist {
                best_dist = d;
                best_pos = Some(pos);
                if d == 0 {
                    break; // Can't do better than exact match
                }
            }
        }

        if best_dist <= max_dist {
            best_pos
        } else {
            None
        }
    }
}

fn find_anchor_candidates(read: &[u8], anchor: &AnchorPlan) -> Vec<AnchorCandidate> {
    let anchor_len = anchor.sequence.len();
    if anchor_len == 0 || anchor_len > read.len() {
        return Vec::new();
    }

    let mut candidates = Vec::new();
    for start in 0..=read.len() - anchor_len {
        let window = &read[start..start + anchor_len];
        let score = match anchor.dist_kind {
            DistanceKind::Hamming => {
                hamming_distance_at_most(window, &anchor.sequence, anchor.max_dist)
            }
            DistanceKind::Levenshtein => u8::MAX,
        };
        if score <= anchor.max_dist {
            candidates.push(AnchorCandidate {
                start,
                end: start + anchor_len,
                score,
            });
        }
    }
    candidates
}

/// Compute Hamming distance between two equal-length byte slices, but stop
/// once the mismatch count exceeds `limit`.
#[inline(always)]
fn hamming_distance_at_most(a: &[u8], b: &[u8], limit: u8) -> u8 {
    debug_assert_eq!(a.len(), b.len());
    let mut dist = 0u8;
    for i in 0..a.len() {
        if a[i] != b[i] {
            dist += 1;
            if dist > limit {
                return dist;
            }
        }
    }
    dist
}

// ── Compilation ─────────────────────────────────────────────────────

/// Pre-scan a FragmentGeom to determine barcode level assignments.
/// Returns a mapping function: given a GeoTagType, returns the barcode level.
fn compute_barcode_levels(geom: &FragmentGeom) -> impl Fn(&GeoTagType) -> Option<u8> {
    let all_parts: Vec<&GeoPart> = geom
        .read1
        .parts
        .iter()
        .chain(geom.read2.parts.iter())
        .collect();

    let has_sample = all_parts
        .iter()
        .any(|p| matches!(p.tag, GeoTagType::SampleBarcode));
    let has_plain_bc = all_parts
        .iter()
        .any(|p| matches!(p.tag, GeoTagType::Barcode));

    // When s + b appear together: s -> level 0, b -> level 1
    // When only numbered barcodes: use their explicit levels
    // When only b (no s): b -> level 0
    let plain_bc_level = if has_sample && has_plain_bc {
        1u8 // b gets level 1 when s is present (s is always level 0)
    } else {
        0u8 // b gets level 0 when s is absent
    };

    move |tag: &GeoTagType| -> Option<u8> {
        match tag {
            GeoTagType::SampleBarcode => Some(0),
            GeoTagType::Barcode => Some(plain_bc_level),
            GeoTagType::NumberedBarcode(n) => Some(*n),
            _ => None,
        }
    }
}

fn compile_boundary_resolved_read_plan(
    read_geom: &ReadGeom,
    bc_level_fn: &dyn Fn(&GeoTagType) -> Option<u8>,
) -> Result<
    (
        BoundaryResolvedReadPlan,
        Vec<usize>,
        Vec<bool>,
        usize,
        bool,
        bool,
    ),
    String,
> {
    let mut anchors = Vec::new();
    let mut segments = Vec::new();
    let mut bc_lens = Vec::new();
    let mut bc_needs_padding = Vec::new();
    let mut umi_len = 0usize;
    let mut umi_needs_padding = false;
    let mut has_read = false;

    let mut current_parts = Vec::new();
    let mut left_boundary = SegmentBoundaryRef::ReadStart;

    for part in &read_geom.parts {
        if matches!(part.tag, GeoTagType::Fixed) {
            let anchor_idx = anchors.len();
            segments.push(ResolvedSegmentPlan {
                left: left_boundary,
                right: SegmentBoundaryRef::AnchorStart(anchor_idx),
                parts: std::mem::take(&mut current_parts),
            });

            anchors.push(AnchorPlan {
                sequence: part.sequence.clone().unwrap_or_default(),
                max_dist: part.tolerance.as_ref().map(|t| t.max_dist).unwrap_or(0),
                dist_kind: part
                    .tolerance
                    .as_ref()
                    .map(|t| t.kind)
                    .unwrap_or(DistanceKind::Hamming),
            });
            left_boundary = SegmentBoundaryRef::AnchorEnd(anchor_idx);
            continue;
        }

        let target = extract_target_for_part(part, bc_level_fn);
        match part.len {
            GeoLen::Fixed(len) => {
                let len = len as usize;
                current_parts.push(ResolvedSegmentPart::Fixed { len, target });
                if let Some(target) = target {
                    record_target_info(
                        target,
                        len,
                        false,
                        &mut bc_lens,
                        &mut bc_needs_padding,
                        &mut umi_len,
                        &mut umi_needs_padding,
                        &mut has_read,
                    );
                }
            }
            GeoLen::Range(min, max) => {
                current_parts.push(ResolvedSegmentPart::Variable {
                    min_len: min as usize,
                    max_len: max as usize,
                    target,
                });
                if let Some(target) = target {
                    record_target_info(
                        target,
                        max as usize,
                        min != max,
                        &mut bc_lens,
                        &mut bc_needs_padding,
                        &mut umi_len,
                        &mut umi_needs_padding,
                        &mut has_read,
                    );
                }
            }
            GeoLen::Unbounded => {
                current_parts.push(ResolvedSegmentPart::Unbounded { target });
                if let Some(target) = target {
                    if matches!(target, ExtractTarget::Read) {
                        has_read = true;
                    }
                }
            }
        }
    }

    segments.push(ResolvedSegmentPlan {
        left: left_boundary,
        right: SegmentBoundaryRef::ReadEnd,
        parts: current_parts,
    });

    Ok((
        BoundaryResolvedReadPlan { anchors, segments },
        bc_lens,
        bc_needs_padding,
        umi_len,
        umi_needs_padding,
        has_read,
    ))
}

/// Compile a ReadGeom into a ReadPlan.
/// Returns (plan, bc_lens_by_level, bc_padding_by_level, umi_len, umi_needs_padding, has_read).
fn compile_read_plan_with_levels(
    read_geom: &ReadGeom,
    bc_level_fn: &dyn Fn(&GeoTagType) -> Option<u8>,
) -> Result<(ReadPlan, Vec<usize>, Vec<bool>, usize, bool, bool), String> {
    let mut steps = Vec::new();
    let mut offset = 0usize;
    let mut bc_lens: Vec<usize> = Vec::new();
    let mut bc_needs_padding: Vec<bool> = Vec::new();
    let mut umi_len = 0usize;
    let mut umi_needs_padding = false;
    let mut has_read = false;
    let mut has_variable = false;
    let mut i = 0;
    while i < read_geom.parts.len() {
        let part = &read_geom.parts[i];

        if matches!(part.len, GeoLen::Range(_, _)) {
            has_variable = true;

            let anchor_idx = read_geom.parts[(i + 1)..]
                .iter()
                .position(|part| matches!(part.tag, GeoTagType::Fixed))
                .map(|rel_idx| i + 1 + rel_idx)
                .ok_or_else(|| {
                    "variable-length fields must be followed by a fixed anchor (f[SEQ])".to_string()
                })?;

            let region_start_offset = offset;
            let mut min_off = offset;
            let mut max_off = offset;
            let mut pre_anchor_parts = Vec::new();
            for pp in &read_geom.parts[i..anchor_idx] {
                let target = extract_target_for_part(pp, bc_level_fn);
                match &pp.len {
                    GeoLen::Fixed(len) => {
                        let len = *len as usize;
                        pre_anchor_parts.push(AnchoredPart::Fixed { len, target });
                        min_off += len;
                        max_off += len;
                        if let Some(target) = target {
                            record_target_len(
                                &mut bc_lens,
                                &mut bc_needs_padding,
                                &mut umi_len,
                                &mut umi_needs_padding,
                                &mut has_read,
                                target,
                                len,
                                false,
                            );
                        }
                    }
                    GeoLen::Range(min, max) => {
                        pre_anchor_parts.push(AnchoredPart::Variable {
                            min_len: *min as usize,
                            max_len: *max as usize,
                            target,
                        });
                        min_off += *min as usize;
                        max_off += *max as usize;
                        if let Some(target) = target {
                            record_target_len(
                                &mut bc_lens,
                                &mut bc_needs_padding,
                                &mut umi_len,
                                &mut umi_needs_padding,
                                &mut has_read,
                                target,
                                *max as usize,
                                min != max,
                            );
                        }
                    }
                    GeoLen::Unbounded => {
                        return Err(
                            "unbounded fields cannot appear before an anchored variable-length region"
                                .into(),
                        );
                    }
                }
            }

            let anchor_part = &read_geom.parts[anchor_idx];
            let anchor_seq = anchor_part.sequence.as_ref().unwrap().clone();
            let (max_dist, dist_kind) = match &anchor_part.tolerance {
                Some(t) => (t.max_dist, t.kind),
                None => (0, DistanceKind::Hamming),
            };

            let mut post_steps = Vec::new();
            let mut post_offset = 0usize;
            for pp in &read_geom.parts[(anchor_idx + 1)..] {
                if matches!(pp.len, GeoLen::Range(_, _)) {
                    return Err(
                        "multiple variable-length anchored regions in the same read are not yet supported"
                            .into(),
                    );
                }
                compile_part_to_step(
                    pp,
                    &mut post_offset,
                    &mut post_steps,
                    &mut bc_lens,
                    &mut bc_needs_padding,
                    &mut umi_len,
                    &mut umi_needs_padding,
                    &mut has_read,
                    bc_level_fn,
                );
            }

            steps.push(ExtractionStep::AnchorSearch {
                region_start_offset,
                min_offset: min_off,
                max_offset: max_off,
                anchor: anchor_seq,
                max_dist,
                dist_kind,
                pre_anchor_parts,
                post_anchor_steps: post_steps,
            });

            break;
        }

        match (&part.tag, &part.len) {
            (GeoTagType::Fixed, GeoLen::Fixed(len)) if part.tolerance.is_none() => {
                offset += *len as usize;
            }
            _ => {
                compile_part_to_step(
                    part,
                    &mut offset,
                    &mut steps,
                    &mut bc_lens,
                    &mut bc_needs_padding,
                    &mut umi_len,
                    &mut umi_needs_padding,
                    &mut has_read,
                    bc_level_fn,
                );
            }
        }

        i += 1;
    }

    Ok((
        ReadPlan {
            steps,
            has_variable,
        },
        bc_lens,
        bc_needs_padding,
        umi_len,
        umi_needs_padding,
        has_read,
    ))
}

/// Compile a single GeoPart into an ExtractionStep.
fn compile_part_to_step(
    part: &GeoPart,
    offset: &mut usize,
    steps: &mut Vec<ExtractionStep>,
    bc_lens: &mut Vec<usize>,
    bc_needs_padding: &mut Vec<bool>,
    umi_len: &mut usize,
    umi_needs_padding: &mut bool,
    has_read: &mut bool,
    bc_level_fn: &dyn Fn(&GeoTagType) -> Option<u8>,
) {
    let target = extract_target_for_part(part, bc_level_fn);

    match &part.len {
        GeoLen::Fixed(len) => {
            let l = *len as usize;
            if let Some(t) = target {
                steps.push(ExtractionStep::FixedSlice {
                    offset: *offset,
                    len: l,
                    target: t,
                });
                record_target_info(
                    t,
                    l,
                    false,
                    bc_lens,
                    bc_needs_padding,
                    umi_len,
                    umi_needs_padding,
                    has_read,
                );
            } else {
                steps.push(ExtractionStep::FixedSkip);
            }
            *offset += l;
        }
        GeoLen::Range(_, max) => {
            // Variable-length barcode/UMI: use max length for the extraction slot.
            // Normalization (padding) happens at a higher level.
            let l = *max as usize;
            if let Some(t) = target {
                steps.push(ExtractionStep::FixedSlice {
                    offset: *offset,
                    len: l,
                    target: t,
                });
                record_target_info(
                    t,
                    l,
                    true,
                    bc_lens,
                    bc_needs_padding,
                    umi_len,
                    umi_needs_padding,
                    has_read,
                );
            }
            // Note: offset advance depends on actual captured length at runtime.
            // For variable-length, this is handled by AnchorSearch in the caller.
            *offset += l;
        }
        GeoLen::Unbounded => {
            if let Some(t) = target {
                steps.push(ExtractionStep::Unbounded {
                    offset: *offset,
                    target: t,
                });
                if matches!(t, ExtractTarget::Read) {
                    *has_read = true;
                }
            } else {
                steps.push(ExtractionStep::UnboundedSkip);
            }
            // Unbounded: no more parts after this
        }
    }
}

fn extract_target_for_part(
    part: &GeoPart,
    bc_level_fn: &dyn Fn(&GeoTagType) -> Option<u8>,
) -> Option<ExtractTarget> {
    match &part.tag {
        GeoTagType::Umi => Some(ExtractTarget::Umi),
        GeoTagType::Read => Some(ExtractTarget::Read),
        GeoTagType::Discard | GeoTagType::Fixed => None,
        tag => bc_level_fn(tag).map(ExtractTarget::Barcode),
    }
}

fn record_target_len(
    bc_lens: &mut Vec<usize>,
    bc_needs_padding: &mut Vec<bool>,
    umi_len: &mut usize,
    umi_needs_padding: &mut bool,
    has_read: &mut bool,
    target: ExtractTarget,
    len: usize,
    needs_padding: bool,
) {
    record_target_info(
        target,
        len,
        needs_padding,
        bc_lens,
        bc_needs_padding,
        umi_len,
        umi_needs_padding,
        has_read,
    );
}

fn record_target_info(
    target: ExtractTarget,
    len: usize,
    needs_padding: bool,
    bc_lens: &mut Vec<usize>,
    bc_needs_padding: &mut Vec<bool>,
    umi_len: &mut usize,
    umi_needs_padding: &mut bool,
    has_read: &mut bool,
) {
    match target {
        ExtractTarget::Barcode(level) => {
            while bc_lens.len() <= level as usize {
                bc_lens.push(0);
                bc_needs_padding.push(false);
            }
            bc_lens[level as usize] = len;
            bc_needs_padding[level as usize] |= needs_padding;
        }
        ExtractTarget::Umi => {
            *umi_len = len;
            *umi_needs_padding |= needs_padding;
        }
        ExtractTarget::Read => *has_read = true,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parse::parse_geometry;

    #[test]
    fn extract_chromium_v3() {
        let geom = parse_geometry("1{b[16]u[12]x:}2{r:}").unwrap();
        let compiled = CompiledGeom::from_fragment_geom(&geom).unwrap();

        let r1 = b"ACGTACGTACGTACGTAAAAAAAAAAAA_extra_data";
        let r2 = b"BIOLOGICAL_READ_SEQUENCE_HERE";

        let result = compiled.extract(r1, r2);
        assert_eq!(result.barcodes[0], Some(&r1[..16]));
        assert_eq!(result.umi, Some(&r1[16..28]));
        assert_eq!(result.reads.len(), 1);
        assert_eq!(result.reads[0], r2.as_slice());
    }

    #[test]
    fn extract_flex_v1() {
        let geom = parse_geometry("1{b[16]u[12]x:}2{r[50]x[18]s[8]x:}").unwrap();
        let compiled = CompiledGeom::from_fragment_geom(&geom).unwrap();

        assert_eq!(compiled.meta().num_bc_levels, 2);

        // Build synthetic reads
        let r1 = b"CELLBARCODEACGTUMI_AAAAAAAA_extra_stuff_here";
        let mut r2 = vec![b'N'; 80];
        // Bio read at 0..50, skip 18, sample BC at 68..76
        r2[68..76].copy_from_slice(b"SAMPLEBC");

        let result = compiled.extract(r1, &r2);
        assert_eq!(result.barcodes.len(), 2);
        assert_eq!(result.barcodes[0], Some(&r2[68..76])); // sample BC (s -> b0)
        assert_eq!(result.barcodes[1], Some(&r1[..16])); // cell BC (b -> b1)
        assert_eq!(result.umi, Some(&r1[16..28]));
        assert_eq!(result.reads[0], &r2[..50]);
    }

    #[test]
    fn extract_flex_v2_exact_anchor() {
        let geom = parse_geometry("1{b[16]u[12]x[0-3]f[TTGCTAGGACCG]s[10]x:}2{r:}").unwrap();
        let compiled = CompiledGeom::from_fragment_geom(&geom).unwrap();

        // Build read with 2bp gap before anchor
        let bc = b"ACGTACGTACGTACGT"; // 16bp
        let umi = b"AAAAAAAAAAAA"; // 12bp
        let gap = b"NN"; // 2bp gap
        let anchor = b"TTGCTAGGACCG"; // 12bp
        let sample = b"SAMPLEBC10"; // 10bp
        let rest = b"extra";

        let mut r1 = Vec::new();
        r1.extend_from_slice(bc);
        r1.extend_from_slice(umi);
        r1.extend_from_slice(gap);
        r1.extend_from_slice(anchor);
        r1.extend_from_slice(sample);
        r1.extend_from_slice(rest);

        let r2 = b"BIOLOGICAL_READ";

        let result = compiled.extract(&r1, r2);
        assert_eq!(result.barcodes[1], Some(bc.as_slice())); // cell BC
        assert_eq!(result.umi, Some(umi.as_slice()));
        assert_eq!(result.barcodes[0], Some(sample.as_slice())); // sample BC
        assert_eq!(result.reads[0], r2.as_slice());
    }

    #[test]
    fn extract_flex_v2_no_gap() {
        let geom = parse_geometry("1{b[16]u[12]x[0-3]f[TTGCTAGGACCG]s[10]x:}2{r:}").unwrap();
        let compiled = CompiledGeom::from_fragment_geom(&geom).unwrap();

        // Build read with 0bp gap (anchor immediately after UMI)
        let bc = b"ACGTACGTACGTACGT";
        let umi = b"AAAAAAAAAAAA";
        let anchor = b"TTGCTAGGACCG";
        let sample = b"SAMPLEBC10";

        let mut r1 = Vec::new();
        r1.extend_from_slice(bc);
        r1.extend_from_slice(umi);
        r1.extend_from_slice(anchor);
        r1.extend_from_slice(sample);

        let r2 = b"BIO";

        let result = compiled.extract(&r1, r2);
        assert_eq!(result.barcodes[0], Some(sample.as_slice()));
    }

    #[test]
    fn extract_flex_v2_anchor_not_found() {
        let geom = parse_geometry("1{b[16]u[12]x[0-3]f[TTGCTAGGACCG]s[10]x:}2{r:}").unwrap();
        let compiled = CompiledGeom::from_fragment_geom(&geom).unwrap();

        // Read with wrong anchor
        let r1 = vec![b'A'; 60];
        let r2 = b"BIO";

        let result = compiled.extract(&r1, r2);
        // Sample BC should be None (anchor not found)
        assert_eq!(result.barcodes[0], None);
        // Cell BC and UMI should still be extracted (before the anchor search)
        assert!(result.barcodes[1].is_some());
        assert!(result.umi.is_some());
    }

    #[test]
    fn extract_hamming_tolerance() {
        let geom =
            parse_geometry("1{b[16]u[12]x[0-3]hamming(f[TTGCTAGGACCG],1)s[10]x:}2{r:}").unwrap();
        let compiled = CompiledGeom::from_fragment_geom(&geom).unwrap();

        // Build read with 1 mismatch in anchor
        let bc = b"ACGTACGTACGTACGT";
        let umi = b"AAAAAAAAAAAA";
        let gap = b"N"; // 1bp gap
        let anchor_mutated = b"TTGCTAGGACCA"; // last G->A (1 mismatch)
        let sample = b"SAMPLEBC10";

        let mut r1 = Vec::new();
        r1.extend_from_slice(bc);
        r1.extend_from_slice(umi);
        r1.extend_from_slice(gap);
        r1.extend_from_slice(anchor_mutated);
        r1.extend_from_slice(sample);

        let r2 = b"BIO";

        let result = compiled.extract(&r1, r2);
        assert_eq!(result.barcodes[0], Some(sample.as_slice()));
    }

    #[test]
    fn extract_variable_barcode_with_anchor() {
        let geom = parse_geometry("1{b[9-10]u[12]f[ACGT]x:}2{r:}").unwrap();
        let compiled = CompiledGeom::from_fragment_geom(&geom).unwrap();
        assert!(compiled.meta().normalization.any_normalization);
        assert_eq!(
            compiled.meta().normalization.bc_needs_padding.as_slice(),
            &[true]
        );
        assert!(!compiled.meta().normalization.umi_needs_padding);

        let bc = b"ACGTACGTA";
        let umi = b"TTTTTTTTTTTT";
        let anchor = b"ACGT";
        let mut r1 = Vec::new();
        r1.extend_from_slice(bc);
        r1.extend_from_slice(umi);
        r1.extend_from_slice(anchor);
        r1.extend_from_slice(b"TAIL");

        let result = compiled.extract(&r1, b"BIO");
        assert_eq!(result.barcodes[0], Some(bc.as_slice()));
        assert_eq!(result.umi, Some(umi.as_slice()));
    }

    #[test]
    fn extract_boundary_resolved_prefix_read_and_suffix_barcode() {
        let geom = parse_geometry("1{r:f[ACAGT]b[9-11]}2{u[12]x:}").unwrap();
        let compiled = CompiledGeom::from_fragment_geom(&geom).unwrap();
        assert!(matches!(compiled, CompiledGeom::BoundaryResolved(_)));

        let read_prefix = b"BIOREAD";
        let anchor = b"ACAGT";
        let barcode = b"BARCODE09";
        let mut r1 = Vec::new();
        r1.extend_from_slice(read_prefix);
        r1.extend_from_slice(anchor);
        r1.extend_from_slice(barcode);

        let r2 = b"TTTTTTTTTTTTtail";
        let result = compiled.extract(&r1, r2);
        assert_eq!(result.reads[0], read_prefix.as_slice());
        assert_eq!(result.barcodes[0], Some(barcode.as_slice()));
        assert_eq!(result.umi, Some(&r2[..12]));
    }

    #[test]
    fn boundary_resolved_prefers_leftmost_best_anchor() {
        let geom = parse_geometry("1{r:f[AC]b[2-7]}2{u[12]x:}").unwrap();
        let compiled = CompiledGeom::from_fragment_geom(&geom).unwrap();

        let r1 = b"READACXYZACBC";
        let r2 = b"TTTTTTTTTTTT";
        let result = compiled.extract(r1, r2);
        assert_eq!(result.reads[0], b"READ".as_slice());
        assert_eq!(result.barcodes[0], Some(b"XYZACBC".as_slice()));
    }

    #[test]
    fn general_extractor_does_not_return_truncated_fields() {
        let geom = parse_geometry("1{b[16]u[12]x[0-3]f[TTGCTAGGACCG]s[10]x:}2{r:}").unwrap();
        let compiled = CompiledGeom::from_fragment_geom(&geom).unwrap();

        let r1 = b"SHORT";
        let r2 = b"BIO";

        let result = compiled.extract(r1, r2);
        assert_eq!(result.barcodes[1], None);
        assert_eq!(result.umi, None);
    }

    #[test]
    fn hamming_distance_test() {
        assert_eq!(hamming_distance_at_most(b"ACGT", b"ACGT", u8::MAX), 0);
        assert_eq!(hamming_distance_at_most(b"ACGT", b"ACGA", u8::MAX), 1);
        assert_eq!(hamming_distance_at_most(b"ACGT", b"TGCA", u8::MAX), 4);
        assert_eq!(hamming_distance_at_most(b"AAAA", b"TTTT", u8::MAX), 4);
    }

    #[test]
    fn hamming_distance_short_circuits_past_limit() {
        assert_eq!(hamming_distance_at_most(b"ACGT", b"TGCA", 1), 2);
        assert_eq!(hamming_distance_at_most(b"ACGT", b"ACGA", 1), 1);
        assert_eq!(hamming_distance_at_most(b"ACGT", b"ACGT", 0), 0);
    }

    #[test]
    fn variable_umi_sets_normalization_metadata() {
        let geom = parse_geometry("1{b[16]u[10-12]f[ACGT]x:}2{r:}").unwrap();
        let compiled = CompiledGeom::from_fragment_geom(&geom).unwrap();
        assert!(compiled.meta().normalization.any_normalization);
        assert_eq!(compiled.meta().bc_lens.as_slice(), &[16]);
        assert_eq!(compiled.meta().umi_len, 12);
        assert_eq!(
            compiled.meta().normalization.bc_needs_padding.as_slice(),
            &[false]
        );
        assert!(compiled.meta().normalization.umi_needs_padding);
    }
}
