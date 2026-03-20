//! Sequence extraction from reads using a compiled geometry.
//!
//! The extraction is split into two phases:
//! 1. **Compile**: Convert a parsed [`FragmentGeom`] into a [`CompiledGeom`] with
//!    precomputed offsets and extraction plans. This happens once at startup.
//! 2. **Extract**: Apply the compiled geometry to each read pair. This is the hot
//!    path and must be as fast as possible — zero allocation, minimal branching.
//!
//! For fixed-length geometries (the common case), extraction is a simple slice
//! operation: `&read[offset..offset+len]`. For variable-length geometries with
//! anchors, a search phase locates the anchor within a small window first.

use smallvec::{smallvec, SmallVec};

use crate::types::*;

/// A precomputed extraction plan for one read.
#[derive(Debug, Clone)]
struct ReadPlan {
    /// Extraction steps, in order. Each step extracts or skips a region.
    steps: Vec<ExtractionStep>,
    /// Whether this read contains any variable-length steps (requires search).
    has_variable: bool,
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
    FixedSkip { offset: usize, len: usize },
    /// Extract from offset to end of read.
    Unbounded {
        offset: usize,
        target: ExtractTarget,
    },
    /// Skip to end of read (unbounded discard).
    UnboundedSkip { offset: usize },
    /// Search for an anchor sequence within a window, then continue
    /// extraction relative to the anchor position.
    AnchorSearch {
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
        /// Steps to execute after the anchor is found. Offsets are relative
        /// to the position immediately after the anchor.
        post_anchor_steps: Vec<ExtractionStep>,
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
#[derive(Debug, Clone)]
pub enum CompiledGeom {
    /// Fast path for simple single-barcode geometries (e.g. Chromium v2/v3).
    /// BC+UMI on one read at fixed offsets, bio read on the other.
    Simple(SimpleExtractor),
    /// General path for multi-barcode or variable-length geometries.
    General(GeneralExtractor),
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
}

/// Fast extractor for simple single-barcode geometries.
#[derive(Debug, Clone)]
pub struct SimpleExtractor {
    pub meta: GeomMeta,
    sg: SimpleGeom,
}

/// General extractor for complex geometries.
#[derive(Debug, Clone)]
pub struct GeneralExtractor {
    pub meta: GeomMeta,
    r1_plan: ReadPlan,
    r2_plan: ReadPlan,
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

impl GeneralExtractor {
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

impl CompiledGeom {
    /// Compile a parsed geometry into an extraction plan.
    pub fn from_fragment_geom(geom: &FragmentGeom) -> Result<Self, String> {
        let bc_level_fn = compute_barcode_levels(geom);
        let (r1_plan, r1_bc_info, r1_umi_len, r1_has_read) =
            compile_read_plan_with_levels(&geom.read1, &bc_level_fn)?;
        let (r2_plan, r2_bc_info, r2_umi_len, r2_has_read) =
            compile_read_plan_with_levels(&geom.read2, &bc_level_fn)?;

        // Merge barcode info from both reads
        let mut bc_lens = SmallVec::new();
        let num_bc_levels = r1_bc_info.len().max(r2_bc_info.len());
        for level in 0..num_bc_levels {
            let len = r1_bc_info
                .get(level)
                .or_else(|| r2_bc_info.get(level))
                .copied()
                .unwrap_or(0);
            bc_lens.push(len);
        }

        let umi_len = r1_umi_len.max(r2_umi_len);

        let is_paired_read = r1_has_read && r2_has_read;

        // Detect simple geometry: 1 barcode level, no variable steps,
        // BC+UMI on R1 at fixed offsets, bio read on R2 (unbounded).
        let simple = if num_bc_levels == 1
            && !r1_plan.has_variable
            && !r2_plan.has_variable
        {
            // Look for BC and UMI fixed slices in R1 plan
            let mut bc_info = None;
            let mut umi_info = None;
            let mut r1_read_offset = None;
            for step in &r1_plan.steps {
                match step {
                    ExtractionStep::FixedSlice { offset, len, target: ExtractTarget::Barcode(0) } => {
                        bc_info = Some((*offset, *len));
                    }
                    ExtractionStep::FixedSlice { offset, len, target: ExtractTarget::Umi } => {
                        umi_info = Some((*offset, *len));
                    }
                    ExtractionStep::Unbounded { offset, target: ExtractTarget::Read } => {
                        r1_read_offset = Some(*offset);
                    }
                    _ => {}
                }
            }
            // Check if R2 has an unbounded read
            let r2_has_unbounded_read = r2_plan.steps.iter().any(|s| {
                matches!(s, ExtractionStep::Unbounded { offset: 0, target: ExtractTarget::Read })
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
        } else {
            None
        };

        let meta = GeomMeta {
            num_bc_levels,
            bc_lens,
            umi_len,
            is_paired_read,
        };

        if let Some(sg) = simple {
            Ok(CompiledGeom::Simple(SimpleExtractor { meta, sg }))
        } else {
            Ok(CompiledGeom::General(GeneralExtractor {
                meta,
                r1_plan,
                r2_plan,
            }))
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
    ///     CompiledGeom::General(ext) => {
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
            CompiledGeom::General(ext) => ext.extract(r1, r2),
        }
    }

    /// Access geometry metadata (barcode lengths, UMI length, etc.).
    pub fn meta(&self) -> &GeomMeta {
        match self {
            CompiledGeom::Simple(ext) => &ext.meta,
            CompiledGeom::General(ext) => &ext.meta,
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

/// Execute a single fixed-offset step (no branching on variable-length).
#[inline(always)]
fn execute_fixed_step<'a>(
    step: &ExtractionStep,
    read: &'a [u8],
    result: &mut ExtractedSeqs<'a>,
) {
    match step {
        ExtractionStep::FixedSlice {
            offset,
            len,
            target,
        } => {
            let end = (*offset + *len).min(read.len());
            if *offset < read.len() {
                let slice = &read[*offset..end];
                assign_target(target, slice, result);
            }
        }
        ExtractionStep::Unbounded { offset, target } => {
            if *offset < read.len() {
                assign_target(target, &read[*offset..], result);
            }
        }
        ExtractionStep::FixedSkip { .. } | ExtractionStep::UnboundedSkip { .. } => {
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
                min_offset,
                max_offset,
                anchor,
                max_dist,
                dist_kind,
                post_anchor_steps,
            } => {
                // Search for the anchor within the window
                if let Some(anchor_pos) =
                    find_anchor(read, *min_offset, *max_offset, anchor, *max_dist, *dist_kind)
                {
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
            let end = (abs_offset + *len).min(read.len());
            if abs_offset < read.len() {
                assign_target(target, &read[abs_offset..end], result);
            }
        }
        ExtractionStep::Unbounded { offset, target } => {
            let abs_offset = base + *offset;
            if abs_offset < read.len() {
                assign_target(target, &read[abs_offset..], result);
            }
        }
        ExtractionStep::FixedSkip { .. } | ExtractionStep::UnboundedSkip { .. } => {}
        ExtractionStep::AnchorSearch { .. } => {
            // Nested anchor searches not supported
        }
    }
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
                    hamming_distance(&read[pos..pos + anchor_len], anchor)
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

/// Compute Hamming distance between two equal-length byte slices.
#[inline(always)]
fn hamming_distance(a: &[u8], b: &[u8]) -> u8 {
    debug_assert_eq!(a.len(), b.len());
    let mut dist = 0u8;
    for i in 0..a.len() {
        if a[i] != b[i] {
            dist += 1;
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

/// Compile a ReadGeom into a ReadPlan.
/// Returns (plan, bc_lens_by_level, umi_len, has_read).
fn compile_read_plan_with_levels(
    read_geom: &ReadGeom,
    bc_level_fn: &dyn Fn(&GeoTagType) -> Option<u8>,
) -> Result<(ReadPlan, Vec<usize>, usize, bool), String> {
    let mut steps = Vec::new();
    let mut offset = 0usize;
    let mut bc_lens: Vec<usize> = Vec::new();
    let mut umi_len = 0usize;
    let mut has_read = false;
    let mut has_variable = false;
    let bc_level_fn = bc_level_fn;

    let mut i = 0;
    while i < read_geom.parts.len() {
        let part = &read_geom.parts[i];

        match (&part.tag, &part.len) {
            // Variable-length discard followed by anchor: compile as AnchorSearch
            (GeoTagType::Discard, GeoLen::Range(min, max)) => {
                has_variable = true;
                let min_off = offset + *min as usize;
                let max_off = offset + *max as usize;

                // The next part must be a Fixed anchor (validated earlier)
                let anchor_part = &read_geom.parts[i + 1];
                let anchor_seq = anchor_part.sequence.as_ref().unwrap().clone();
                let (max_dist, dist_kind) = match &anchor_part.tolerance {
                    Some(t) => (t.max_dist, t.kind),
                    None => (0, DistanceKind::Hamming),
                };
                let anchor_len = anchor_seq.len();

                // Compile post-anchor steps with relative offsets (from end of anchor)
                let mut post_steps = Vec::new();
                let mut post_offset = 0usize;
                for j in (i + 2)..read_geom.parts.len() {
                    let pp = &read_geom.parts[j];
                    compile_part_to_step(
                        pp,
                        &mut post_offset,
                        &mut post_steps,
                        &mut bc_lens,
                        &mut umi_len,
                        &mut has_read,
                        bc_level_fn,
                    );
                }

                steps.push(ExtractionStep::AnchorSearch {
                    min_offset: min_off,
                    max_offset: max_off,
                    anchor: anchor_seq,
                    max_dist,
                    dist_kind,
                    post_anchor_steps: post_steps,
                });

                // Skip all remaining parts (handled by post_anchor_steps)
                break;
            }

            // Fixed anchor (not preceded by variable discard): just skip it
            (GeoTagType::Fixed, GeoLen::Fixed(len)) if part.tolerance.is_none() => {
                offset += *len as usize;
            }

            // All other parts: compile to fixed steps
            _ => {
                compile_part_to_step(
                    part,
                    &mut offset,
                    &mut steps,
                    &mut bc_lens,
                    &mut umi_len,
                    &mut has_read,
                    bc_level_fn,
                );
            }
        }

        i += 1;
    }

    Ok((ReadPlan { steps, has_variable }, bc_lens, umi_len, has_read))
}

/// Compile a single GeoPart into an ExtractionStep.
fn compile_part_to_step(
    part: &GeoPart,
    offset: &mut usize,
    steps: &mut Vec<ExtractionStep>,
    bc_lens: &mut Vec<usize>,
    umi_len: &mut usize,
    has_read: &mut bool,
    bc_level_fn: &dyn Fn(&GeoTagType) -> Option<u8>,
) {
    let target = match &part.tag {
        GeoTagType::Umi => Some(ExtractTarget::Umi),
        GeoTagType::Read => Some(ExtractTarget::Read),
        GeoTagType::Discard | GeoTagType::Fixed => None,
        tag => bc_level_fn(tag).map(ExtractTarget::Barcode),
    };

    match &part.len {
        GeoLen::Fixed(len) => {
            let l = *len as usize;
            if let Some(t) = target {
                steps.push(ExtractionStep::FixedSlice {
                    offset: *offset,
                    len: l,
                    target: t,
                });
                match t {
                    ExtractTarget::Barcode(level) => {
                        while bc_lens.len() <= level as usize {
                            bc_lens.push(0);
                        }
                        bc_lens[level as usize] = l;
                    }
                    ExtractTarget::Umi => *umi_len = l,
                    ExtractTarget::Read => *has_read = true,
                }
            } else {
                steps.push(ExtractionStep::FixedSkip {
                    offset: *offset,
                    len: l,
                });
            }
            *offset += l;
        }
        GeoLen::Range(min, max) => {
            // Variable-length barcode/UMI: use max length for the extraction slot.
            // Normalization (padding) happens at a higher level.
            let l = *max as usize;
            if let Some(t) = target {
                steps.push(ExtractionStep::FixedSlice {
                    offset: *offset,
                    len: l,
                    target: t,
                });
                match t {
                    ExtractTarget::Barcode(level) => {
                        while bc_lens.len() <= level as usize {
                            bc_lens.push(0);
                        }
                        bc_lens[level as usize] = l;
                    }
                    ExtractTarget::Umi => *umi_len = l,
                    ExtractTarget::Read => *has_read = true,
                }
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
                steps.push(ExtractionStep::UnboundedSkip { offset: *offset });
            }
            // Unbounded: no more parts after this
        }
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
        let geom =
            parse_geometry("1{b[16]u[12]x[0-3]f[TTGCTAGGACCG]s[10]x:}2{r:}").unwrap();
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
        let geom =
            parse_geometry("1{b[16]u[12]x[0-3]f[TTGCTAGGACCG]s[10]x:}2{r:}").unwrap();
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
        let geom =
            parse_geometry("1{b[16]u[12]x[0-3]f[TTGCTAGGACCG]s[10]x:}2{r:}").unwrap();
        let compiled = CompiledGeom::from_fragment_geom(&geom).unwrap();

        // Read with wrong anchor
        let mut r1 = vec![b'A'; 60];
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
        let geom = parse_geometry(
            "1{b[16]u[12]x[0-3]hamming(f[TTGCTAGGACCG],1)s[10]x:}2{r:}",
        )
        .unwrap();
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
    fn hamming_distance_test() {
        assert_eq!(hamming_distance(b"ACGT", b"ACGT"), 0);
        assert_eq!(hamming_distance(b"ACGT", b"ACGA"), 1);
        assert_eq!(hamming_distance(b"ACGT", b"TGCA"), 4);
        assert_eq!(hamming_distance(b"AAAA", b"TTTT"), 4);
    }
}
