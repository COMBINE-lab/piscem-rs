//! Chumsky-based parser for sequencing read geometry descriptions.
//!
//! Grammar:
//! ```text
//! FragmentGeom := Read1Desc Read2Desc
//! Read1Desc    := '1{' DescList '}'
//! Read2Desc    := '2{' DescList '}'
//! DescList     := Desc+
//! Desc         := TagSpec | DistFunc | UnboundedDesc
//! TagSpec      := TagChar '[' LenSpec ']'
//!              |  TagChar '[' Sequence ']'    (* for 'f' tag *)
//!              |  TagChar Digit '[' LenSpec ']'  (* numbered barcode *)
//! LenSpec      := Number ('-' Number)?
//! UnboundedDesc := ('x' | 'r') ':'
//! DistFunc     := 'hamming(' TagSpec ',' Number ')'
//! Number       := [1-9][0-9]* | '0'
//! Sequence     := [ACGT]+
//! TagChar      := 'b' | 'u' | 'r' | 'f' | 'x' | 's'
//! ```

use chumsky::prelude::*;

use crate::types::*;

/// Parse a geometry description string into a [`FragmentGeom`].
///
/// Returns `Ok(FragmentGeom)` on success, or an error with source spans
/// and descriptive messages on failure.
pub fn parse_geometry(input: &str) -> Result<FragmentGeom, Vec<Rich<'_, char>>> {
    let result = fragment_geom_parser().parse(input);
    result.into_result()
}

/// Format parse errors into human-readable strings.
pub fn format_errors(input: &str, errors: &[Rich<'_, char>]) -> String {
    use ariadne::{Color, Label, Report, ReportKind, Source};

    let mut buf = Vec::new();
    for err in errors {
        let span = err.span();
        Report::build(ReportKind::Error, span.start..span.end)
            .with_message(format!("{}", err))
            .with_label(
                Label::new(span.start..span.end)
                    .with_message(format!("{}", err.reason()))
                    .with_color(Color::Red),
            )
            .finish()
            .write(Source::from(input), &mut buf)
            .unwrap();
    }
    String::from_utf8(buf).unwrap_or_else(|_| "error rendering failed".to_string())
}

// ── Parser combinators ──────────────────────────────────────────────

fn fragment_geom_parser<'a>() -> impl Parser<'a, &'a str, FragmentGeom, extra::Err<Rich<'a, char>>>
{
    read_desc_parser('1')
        .then(read_desc_parser('2'))
        .then_ignore(end())
        .map(|(r1, r2)| FragmentGeom {
            read1: r1,
            read2: r2,
        })
}

fn read_desc_parser<'a>(
    read_num: char,
) -> impl Parser<'a, &'a str, ReadGeom, extra::Err<Rich<'a, char>>> {
    just(read_num)
        .ignore_then(desc_list_parser().delimited_by(just('{'), just('}')))
        .map(|parts| ReadGeom { parts })
}

fn desc_list_parser<'a>() -> impl Parser<'a, &'a str, Vec<GeoPart>, extra::Err<Rich<'a, char>>> {
    desc_parser().repeated().at_least(1).collect::<Vec<_>>()
}

fn desc_parser<'a>() -> impl Parser<'a, &'a str, GeoPart, extra::Err<Rich<'a, char>>> {
    dist_func_parser()
        .or(tag_spec_parser())
        .or(unbounded_desc_parser())
}

fn tag_spec_parser<'a>() -> impl Parser<'a, &'a str, GeoPart, extra::Err<Rich<'a, char>>> {
    // Numbered barcode: b0[N], b1[N], etc.
    let numbered_bc = just('b')
        .ignore_then(one_of("0123456789").map(|c: char| c.to_digit(10).unwrap() as u8))
        .then(len_spec_parser().delimited_by(just('['), just(']')))
        .map(|(level, len)| GeoPart {
            tag: GeoTagType::NumberedBarcode(level),
            len,
            sequence: None,
            tolerance: None,
        });

    // Fixed sequence: f[ACGT...]
    let fixed_seq = just('f')
        .ignore_then(
            one_of("ACGTacgt")
                .repeated()
                .at_least(1)
                .collect::<String>()
                .delimited_by(just('['), just(']')),
        )
        .map(|seq| GeoPart {
            tag: GeoTagType::Fixed,
            len: GeoLen::Fixed(seq.len() as u32),
            sequence: Some(seq.to_uppercase().into_bytes()),
            tolerance: None,
        });

    // Standard tags: b[N], u[N], r[N], x[N], x[N-M], s[N]
    let standard_tag = one_of("buxrs")
        .then(len_spec_parser().delimited_by(just('['), just(']')))
        .map(|(tag_char, len)| {
            let tag = match tag_char {
                'b' => GeoTagType::Barcode,
                'u' => GeoTagType::Umi,
                'x' => GeoTagType::Discard,
                'r' => GeoTagType::Read,
                's' => GeoTagType::SampleBarcode,
                _ => unreachable!(),
            };
            GeoPart {
                tag,
                len,
                sequence: None,
                tolerance: None,
            }
        });

    // Try numbered barcode first (b0[...]), then fixed, then standard
    numbered_bc.or(fixed_seq).or(standard_tag)
}

fn len_spec_parser<'a>() -> impl Parser<'a, &'a str, GeoLen, extra::Err<Rich<'a, char>>> {
    let num = text::int(10).map(|s: &str| s.parse::<u32>().unwrap());

    num.then(just('-').ignore_then(num).or_not())
        .map(|(n, m)| match m {
            Some(m) => GeoLen::Range(n, m),
            None => GeoLen::Fixed(n),
        })
}

fn unbounded_desc_parser<'a>() -> impl Parser<'a, &'a str, GeoPart, extra::Err<Rich<'a, char>>> {
    one_of("xr").then_ignore(just(':')).map(|tag_char| {
        let tag = match tag_char {
            'x' => GeoTagType::Discard,
            'r' => GeoTagType::Read,
            _ => unreachable!(),
        };
        GeoPart {
            tag,
            len: GeoLen::Unbounded,
            sequence: None,
            tolerance: None,
        }
    })
}

fn dist_func_parser<'a>() -> impl Parser<'a, &'a str, GeoPart, extra::Err<Rich<'a, char>>> {
    let num = text::int(10).map(|s: &str| s.parse::<u8>().unwrap());

    // hamming(f[SEQ], N) or hamming(tag_spec, N)
    let hamming = just("hamming")
        .ignore_then(
            just('(')
                .ignore_then(tag_spec_parser())
                .then_ignore(just(','))
                .then_ignore(just(' ').or_not())
                .then(num)
                .then_ignore(just(')')),
        )
        .map(|(mut part, max_dist)| {
            part.tolerance = Some(MatchTolerance {
                kind: DistanceKind::Hamming,
                max_dist,
            });
            part
        });

    // Future: levenshtein(tag_spec, N)
    // let levenshtein = ...

    hamming
}

// ── Validation ──────────────────────────────────────────────────────

/// Validate a parsed geometry for semantic correctness.
pub fn validate_geometry(geom: &FragmentGeom) -> Result<(), String> {
    let all_parts: Vec<&GeoPart> = geom
        .read1
        .parts
        .iter()
        .chain(geom.read2.parts.iter())
        .collect();

    // Must have at least one barcode
    let has_barcode = all_parts.iter().any(|p| {
        matches!(
            p.tag,
            GeoTagType::Barcode | GeoTagType::SampleBarcode | GeoTagType::NumberedBarcode(_)
        )
    });
    if !has_barcode {
        return Err("geometry must contain at least one barcode tag (b, s, or bN)".into());
    }

    // Must have at least one UMI
    let has_umi = all_parts.iter().any(|p| matches!(p.tag, GeoTagType::Umi));
    if !has_umi {
        return Err("geometry must contain at least one UMI tag (u)".into());
    }

    // Must have at least one read
    let has_read = all_parts.iter().any(|p| matches!(p.tag, GeoTagType::Read));
    if !has_read {
        return Err("geometry must contain at least one biological read tag (r)".into());
    }

    // Validate range widths
    for part in &all_parts {
        if let GeoLen::Range(min, max) = &part.len {
            if min > max {
                return Err(format!(
                    "range minimum {} exceeds maximum {} in {:?} tag",
                    min, max, part.tag
                ));
            }
            if max - min > MAX_RANGE_WIDTH {
                return Err(format!(
                    "range width {} exceeds maximum {} in {:?} tag",
                    max - min,
                    MAX_RANGE_WIDTH,
                    part.tag
                ));
            }
        }
    }

    for read in [&geom.read1, &geom.read2] {
        match classify_read_complexity(read) {
            GeometryComplexity::FixedOffsets => {}
            GeometryComplexity::InferableVariable => validate_inferable_read(read)?,
            GeometryComplexity::BoundaryResolved => validate_boundary_resolved_read(read)?,
        }
    }

    Ok(())
}

/// Classify the executor complexity tier required for a geometry.
///
/// The classifier chooses the most complex tier required by either read:
/// - [`GeometryComplexity::FixedOffsets`] for fully static layouts such as
///   `1{b[16]u[12]x:}2{r:}`
/// - [`GeometryComplexity::InferableVariable`] for a single right-bounded
///   variable-width region such as `1{b[9-10]f[ACGT]u[12]}2{r:}`
/// - [`GeometryComplexity::BoundaryResolved`] for geometries with interior
///   unbounded segments such as `1{r:f[ACAGT]b[9-11]}2{u[12]x:}`
///
/// This is the same tier selection used by [`crate::extract::CompiledGeom`].
pub fn geometry_complexity(geom: &FragmentGeom) -> GeometryComplexity {
    classify_read_complexity(&geom.read1).max(classify_read_complexity(&geom.read2))
}

fn classify_read_complexity(read: &ReadGeom) -> GeometryComplexity {
    let has_non_terminal_unbounded = read
        .parts
        .iter()
        .enumerate()
        .any(|(i, part)| matches!(part.len, GeoLen::Unbounded) && i + 1 < read.parts.len());
    if has_non_terminal_unbounded {
        return GeometryComplexity::BoundaryResolved;
    }

    let has_variable = read
        .parts
        .iter()
        .any(|part| matches!(part.len, GeoLen::Range(_, _)));
    if has_variable {
        GeometryComplexity::InferableVariable
    } else {
        GeometryComplexity::FixedOffsets
    }
}

fn validate_inferable_read(read: &ReadGeom) -> Result<(), String> {
    let Some(first_range_idx) = read
        .parts
        .iter()
        .position(|part| matches!(part.len, GeoLen::Range(_, _)))
    else {
        return Ok(());
    };

    let Some(anchor_rel_idx) = read.parts[(first_range_idx + 1)..]
        .iter()
        .position(|part| matches!(part.tag, GeoTagType::Fixed))
    else {
        return Err(
            "variable-length fields must be followed by a fixed anchor (f[SEQ]) so their boundaries can be inferred"
                .into(),
        );
    };

    let anchor_idx = first_range_idx + 1 + anchor_rel_idx;
    let variable_count = read.parts[first_range_idx..anchor_idx]
        .iter()
        .filter(|part| matches!(part.len, GeoLen::Range(_, _)))
        .count();
    if variable_count > 1 {
        return Err(
            "at most one variable-length field before a fixed anchor is currently supported".into(),
        );
    }

    if read.parts[(anchor_idx + 1)..]
        .iter()
        .any(|part| matches!(part.len, GeoLen::Range(_, _)))
    {
        return Err(
            "multiple variable-length anchored regions in the same read are not yet supported"
                .into(),
        );
    }

    Ok(())
}

fn validate_boundary_resolved_read(read: &ReadGeom) -> Result<(), String> {
    let mut segment_start = 0usize;
    let mut saw_anchor = false;
    for (idx, part) in read.parts.iter().enumerate() {
        if matches!(part.tag, GeoTagType::Fixed) {
            validate_boundary_segment(&read.parts[segment_start..idx])?;
            segment_start = idx + 1;
            saw_anchor = true;
        }
    }

    validate_boundary_segment(&read.parts[segment_start..])?;

    let has_non_terminal_unbounded = read
        .parts
        .iter()
        .enumerate()
        .any(|(i, part)| matches!(part.len, GeoLen::Unbounded) && i + 1 < read.parts.len());
    if has_non_terminal_unbounded && !saw_anchor {
        return Err(
            "non-terminal unbounded fields require at least one fixed anchor (f[SEQ]) to resolve boundaries"
                .into(),
        );
    }

    Ok(())
}

fn validate_boundary_segment(parts: &[GeoPart]) -> Result<(), String> {
    let flexible_count = parts
        .iter()
        .filter(|part| matches!(part.len, GeoLen::Range(_, _) | GeoLen::Unbounded))
        .count();
    if flexible_count > 1 {
        return Err(
            "at most one variable-width or unbounded field is supported within each boundary-resolved segment"
                .into(),
        );
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_chromium_v3() {
        let geom = parse_geometry("1{b[16]u[12]x:}2{r:}").unwrap();
        assert_eq!(geom.read1.parts.len(), 3);
        assert_eq!(geom.read1.parts[0].tag, GeoTagType::Barcode);
        assert_eq!(geom.read1.parts[0].len, GeoLen::Fixed(16));
        assert_eq!(geom.read1.parts[1].tag, GeoTagType::Umi);
        assert_eq!(geom.read1.parts[1].len, GeoLen::Fixed(12));
        assert_eq!(geom.read1.parts[2].tag, GeoTagType::Discard);
        assert_eq!(geom.read1.parts[2].len, GeoLen::Unbounded);
        assert_eq!(geom.read2.parts.len(), 1);
        assert_eq!(geom.read2.parts[0].tag, GeoTagType::Read);
        assert_eq!(geom.read2.parts[0].len, GeoLen::Unbounded);
    }

    #[test]
    fn parse_flex_v1() {
        let geom = parse_geometry("1{b[16]u[12]x:}2{r[50]x[18]s[8]x:}").unwrap();
        assert_eq!(geom.read2.parts.len(), 4);
        assert_eq!(geom.read2.parts[0].tag, GeoTagType::Read);
        assert_eq!(geom.read2.parts[0].len, GeoLen::Fixed(50));
        assert_eq!(geom.read2.parts[2].tag, GeoTagType::SampleBarcode);
        assert_eq!(geom.read2.parts[2].len, GeoLen::Fixed(8));
    }

    #[test]
    fn parse_flex_v2_with_anchor() {
        let geom = parse_geometry("1{b[16]u[12]x[0-3]f[TTGCTAGGACCG]s[10]x:}2{r:}").unwrap();
        assert_eq!(geom.read1.parts.len(), 6);
        assert_eq!(geom.read1.parts[2].tag, GeoTagType::Discard);
        assert_eq!(geom.read1.parts[2].len, GeoLen::Range(0, 3));
        assert_eq!(geom.read1.parts[3].tag, GeoTagType::Fixed);
        assert_eq!(geom.read1.parts[3].sequence, Some(b"TTGCTAGGACCG".to_vec()));
        assert_eq!(geom.read1.parts[4].tag, GeoTagType::SampleBarcode);
        assert_eq!(geom.read1.parts[4].len, GeoLen::Fixed(10));
    }

    #[test]
    fn parse_hamming_distance() {
        let geom =
            parse_geometry("1{b[16]u[12]x[0-3]hamming(f[TTGCTAGGACCG],1)s[10]x:}2{r:}").unwrap();
        let anchor = &geom.read1.parts[3];
        assert_eq!(anchor.tag, GeoTagType::Fixed);
        assert_eq!(
            anchor.tolerance,
            Some(MatchTolerance {
                kind: DistanceKind::Hamming,
                max_dist: 1,
            })
        );
    }

    #[test]
    fn parse_numbered_barcodes() {
        let geom = parse_geometry("1{b0[8]b1[16]u[12]x:}2{r:}").unwrap();
        assert_eq!(geom.read1.parts[0].tag, GeoTagType::NumberedBarcode(0));
        assert_eq!(geom.read1.parts[0].len, GeoLen::Fixed(8));
        assert_eq!(geom.read1.parts[1].tag, GeoTagType::NumberedBarcode(1));
        assert_eq!(geom.read1.parts[1].len, GeoLen::Fixed(16));
    }

    #[test]
    fn parse_fixed_sequence() {
        let geom = parse_geometry("1{b[16]u[12]f[ACGT]r:}2{r:}").unwrap();
        assert_eq!(geom.read1.parts[2].tag, GeoTagType::Fixed);
        assert_eq!(geom.read1.parts[2].sequence, Some(b"ACGT".to_vec()));
        assert_eq!(geom.read1.parts[2].len, GeoLen::Fixed(4));
    }

    #[test]
    fn parse_boundary_resolved_unbounded() {
        let geom = parse_geometry("1{r:f[ACAGT]b[9-11]}2{u[12]x:}").unwrap();
        assert_eq!(geom.read1.parts.len(), 3);
        assert_eq!(geom.read1.parts[0].len, GeoLen::Unbounded);
        assert_eq!(geom.read1.parts[1].tag, GeoTagType::Fixed);
        assert_eq!(geom.read1.parts[2].len, GeoLen::Range(9, 11));
    }

    #[test]
    fn validate_missing_barcode() {
        let geom = parse_geometry("1{u[12]x:}2{r:}").unwrap();
        assert!(validate_geometry(&geom).is_err());
    }

    #[test]
    fn validate_missing_umi() {
        let geom = parse_geometry("1{b[16]x:}2{r:}").unwrap();
        assert!(validate_geometry(&geom).is_err());
    }

    #[test]
    fn validate_range_without_anchor() {
        let geom = parse_geometry("1{b[16]u[12]x[0-3]s[10]x:}2{r:}").unwrap();
        let err = validate_geometry(&geom);
        assert!(err.is_err());
        assert!(err
            .unwrap_err()
            .contains("must be followed by a fixed anchor"));
    }

    #[test]
    fn validate_range_too_wide() {
        let geom = parse_geometry("1{b[16]u[12]x[0-5]f[TTGCTAGGACCG]s[10]x:}2{r:}").unwrap();
        let err = validate_geometry(&geom);
        assert!(err.is_err());
        assert!(err.unwrap_err().contains("exceeds maximum"));
    }

    #[test]
    fn validate_variable_length_barcode_with_anchor() {
        let geom = parse_geometry("1{b[9-10]u[12]f[ACGT]x:}2{r:}").unwrap();
        assert!(validate_geometry(&geom).is_ok());
    }

    #[test]
    fn validate_variable_length_barcode_without_anchor_rejected() {
        let geom = parse_geometry("1{b[9-10]u[12]x:}2{r:}").unwrap();
        let err = validate_geometry(&geom);
        assert!(err.is_err());
        assert!(err
            .unwrap_err()
            .contains("must be followed by a fixed anchor"));
    }

    #[test]
    fn validate_boundary_segment_with_two_flexible_fields_rejected() {
        let geom = parse_geometry("1{r:b[9-11]f[ACAGT]}2{u[12]x:}").unwrap();
        let err = validate_geometry(&geom);
        assert!(err.is_err());
        assert!(err
            .unwrap_err()
            .contains("at most one variable-width or unbounded field"));
    }

    #[test]
    fn error_message_quality() {
        let result = parse_geometry("1{b[16]u[12]z:}2{r:}");
        assert!(result.is_err());
        let errors = result.unwrap_err();
        let msg = format_errors("1{b[16]u[12]z:}2{r:}", &errors);
        // The error message should contain some indication of the problem
        assert!(!msg.is_empty());
        // It should reference the error position (the 'z' character)
        assert!(msg.contains("Error"));
    }

    #[test]
    fn parse_hamming_with_space() {
        // Allow optional space after comma
        let geom =
            parse_geometry("1{b[16]u[12]x[0-3]hamming(f[TTGCTAGGACCG], 1)s[10]x:}2{r:}").unwrap();
        let anchor = &geom.read1.parts[3];
        assert_eq!(
            anchor.tolerance,
            Some(MatchTolerance {
                kind: DistanceKind::Hamming,
                max_dist: 1,
            })
        );
    }

    #[test]
    fn classify_fixed_vs_variable_complexity() {
        let fixed = parse_geometry("1{b[16]u[12]x:}2{r:}").unwrap();
        assert_eq!(
            geometry_complexity(&fixed),
            GeometryComplexity::FixedOffsets
        );

        let variable = parse_geometry("1{b[9-10]u[12]f[ACGT]x:}2{r:}").unwrap();
        assert_eq!(
            geometry_complexity(&variable),
            GeometryComplexity::InferableVariable
        );

        let boundary = parse_geometry("1{r:f[ACAGT]b[9-11]}2{u[12]x:}").unwrap();
        assert_eq!(
            geometry_complexity(&boundary),
            GeometryComplexity::BoundaryResolved
        );
    }
}
