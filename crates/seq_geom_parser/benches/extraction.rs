//! Performance benchmarks for geometry extraction.
//!
//! These benchmarks ensure extraction stays fast (the hot path runs
//! hundreds of millions to billions of times per dataset). Run with:
//!
//!   cargo bench -p seq_geom_parser
//!
//! Results are saved to `target/criterion/` with HTML reports.

use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId, Throughput};
use seq_geom_parser::{parse_geometry, validate_geometry, CompiledGeom, SimpleExtractor, ExtractionBuf, Extractor};

/// Generate a synthetic read of the given length filled with random-ish bases.
fn make_read(len: usize, seed: u8) -> Vec<u8> {
    let bases = [b'A', b'C', b'G', b'T'];
    (0..len).map(|i| bases[((i as u8).wrapping_add(seed)) as usize % 4]).collect()
}

/// Build a Flex v1 R2 read: 50bp bio + 18bp linker + 8bp sample BC + rest
fn make_flex_v1_r2(sample_bc: &[u8; 8]) -> Vec<u8> {
    let mut r2 = make_read(50, 42); // bio read
    r2.extend_from_slice(&[b'N'; 18]); // linker
    r2.extend_from_slice(sample_bc);
    r2.extend_from_slice(&[b'N'; 10]); // extra
    r2
}

/// Build a Flex v2 R1 read: 16bp BC + 12bp UMI + gap + anchor + 10bp sample BC + rest
fn make_flex_v2_r1(gap_len: usize, anchor: &[u8], sample_bc: &[u8; 10]) -> Vec<u8> {
    let mut r1 = make_read(16, 1); // cell BC
    r1.extend(make_read(12, 2)); // UMI
    r1.extend(vec![b'N'; gap_len]); // variable gap
    r1.extend_from_slice(anchor); // anchor
    r1.extend_from_slice(sample_bc);
    r1.extend(vec![b'N'; 20]); // extra
    r1
}

fn bench_parse(c: &mut Criterion) {
    let mut group = c.benchmark_group("parse");

    let geometries = [
        ("chromium_v3", "1{b[16]u[12]x:}2{r:}"),
        ("flex_v1", "1{b[16]u[12]x:}2{r[50]x[18]s[8]x:}"),
        ("flex_v2", "1{b[16]u[12]x[0-3]f[TTGCTAGGACCG]s[10]x:}2{r:}"),
        ("flex_v2_hamming", "1{b[16]u[12]x[0-3]hamming(f[TTGCTAGGACCG],1)s[10]x:}2{r:}"),
    ];

    for (name, geom_str) in &geometries {
        group.bench_with_input(BenchmarkId::new("parse", name), geom_str, |b, g| {
            b.iter(|| parse_geometry(black_box(g)).unwrap());
        });
    }

    group.finish();
}

fn bench_compile(c: &mut Criterion) {
    let mut group = c.benchmark_group("compile");

    let geometries = [
        ("chromium_v3", "1{b[16]u[12]x:}2{r:}"),
        ("flex_v1", "1{b[16]u[12]x:}2{r[50]x[18]s[8]x:}"),
        ("flex_v2", "1{b[16]u[12]x[0-3]f[TTGCTAGGACCG]s[10]x:}2{r:}"),
    ];

    for (name, geom_str) in &geometries {
        let geom = parse_geometry(geom_str).unwrap();
        group.bench_with_input(BenchmarkId::new("compile", name), &geom, |b, g| {
            b.iter(|| CompiledGeom::from_fragment_geom(black_box(g)).unwrap());
        });
    }

    group.finish();
}

fn bench_extract_fixed(c: &mut Criterion) {
    let mut group = c.benchmark_group("extract_fixed");
    group.throughput(Throughput::Elements(1));

    // Chromium v3: simplest case (all fixed offsets)
    let geom = parse_geometry("1{b[16]u[12]x:}2{r:}").unwrap();
    let compiled = CompiledGeom::from_fragment_geom(&geom).unwrap();
    let r1 = make_read(150, 1);
    let r2 = make_read(150, 2);

    group.bench_function("chromium_v3", |b| {
        b.iter(|| compiled.extract(black_box(&r1), black_box(&r2)));
    });

    // Flex v1: two reads with fixed offsets
    let geom = parse_geometry("1{b[16]u[12]x:}2{r[50]x[18]s[8]x:}").unwrap();
    let compiled = CompiledGeom::from_fragment_geom(&geom).unwrap();
    let r1 = make_read(150, 1);
    let r2 = make_flex_v1_r2(b"ACTTTAGG");

    group.bench_function("flex_v1", |b| {
        b.iter(|| compiled.extract(black_box(&r1), black_box(&r2)));
    });

    group.finish();
}

fn bench_extract_anchor(c: &mut Criterion) {
    let mut group = c.benchmark_group("extract_anchor");
    group.throughput(Throughput::Elements(1));

    let anchor = b"TTGCTAGGACCG";
    let sample_bc = b"SAMPLEBC10";

    // Flex v2 with exact anchor match
    let geom = parse_geometry("1{b[16]u[12]x[0-3]f[TTGCTAGGACCG]s[10]x:}2{r:}").unwrap();
    let compiled = CompiledGeom::from_fragment_geom(&geom).unwrap();
    let r2 = make_read(150, 2);

    for gap in [0, 1, 2, 3] {
        let r1 = make_flex_v2_r1(gap, anchor, sample_bc);
        group.bench_with_input(
            BenchmarkId::new("exact_gap", gap),
            &(r1, r2.clone()),
            |b, (r1, r2)| {
                b.iter(|| compiled.extract(black_box(r1), black_box(r2)));
            },
        );
    }

    // Flex v2 with hamming(1) tolerance
    let geom = parse_geometry(
        "1{b[16]u[12]x[0-3]hamming(f[TTGCTAGGACCG],1)s[10]x:}2{r:}",
    )
    .unwrap();
    let compiled_h1 = CompiledGeom::from_fragment_geom(&geom).unwrap();

    // Exact match with hamming enabled (should still be fast)
    let r1_exact = make_flex_v2_r1(2, anchor, sample_bc);
    group.bench_function("hamming1_exact", |b| {
        b.iter(|| compiled_h1.extract(black_box(&r1_exact), black_box(&r2)));
    });

    // 1 mismatch in anchor
    let mut anchor_mut = anchor.to_vec();
    anchor_mut[11] = b'A'; // last base mutated
    let r1_mismatch = make_flex_v2_r1(2, &anchor_mut, sample_bc);
    group.bench_function("hamming1_mismatch", |b| {
        b.iter(|| compiled_h1.extract(black_box(&r1_mismatch), black_box(&r2)));
    });

    // Anchor not found (worst case: searches all positions)
    let r1_nofind = make_read(80, 99);
    group.bench_function("anchor_not_found", |b| {
        b.iter(|| compiled.extract(black_box(&r1_nofind), black_box(&r2)));
    });

    group.finish();
}

/// Hardcoded extraction mimicking ChromiumProtocol::extract_tech_seqs exactly.
/// This is the absolute minimum work: two bounds checks, two slices, return.
#[inline(never)]
fn hardcoded_chromium_v3_extract<'a>(r1: &'a [u8], r2: &'a [u8]) -> (Option<&'a [u8]>, Option<&'a [u8]>, &'a [u8]) {
    let bc_len = 16;
    let umi_len = 12;
    let barcode = if r1.len() >= bc_len { Some(&r1[..bc_len]) } else { None };
    let umi = if r1.len() >= bc_len + umi_len { Some(&r1[bc_len..bc_len + umi_len]) } else { None };
    (barcode, umi, r2)
}

/// Same but for Chromium v2 (16bp BC, 10bp UMI).
#[inline(never)]
fn hardcoded_chromium_v2_extract<'a>(r1: &'a [u8], r2: &'a [u8]) -> (Option<&'a [u8]>, Option<&'a [u8]>, &'a [u8]) {
    let bc_len = 16;
    let umi_len = 10;
    let barcode = if r1.len() >= bc_len { Some(&r1[..bc_len]) } else { None };
    let umi = if r1.len() >= bc_len + umi_len { Some(&r1[bc_len..bc_len + umi_len]) } else { None };
    (barcode, umi, r2)
}

fn bench_hardcoded_vs_compiled(c: &mut Criterion) {
    let mut group = c.benchmark_group("hardcoded_vs_compiled");
    group.throughput(Throughput::Elements(1));

    let r1 = make_read(150, 1);
    let r2 = make_read(150, 2);

    // Hardcoded v3
    group.bench_function("hardcoded_v3", |b| {
        b.iter(|| hardcoded_chromium_v3_extract(black_box(&r1), black_box(&r2)));
    });

    // Compiled v3 via enum dispatch
    let geom_v3 = parse_geometry("1{b[16]u[12]x:}2{r:}").unwrap();
    let compiled_v3 = CompiledGeom::from_fragment_geom(&geom_v3).unwrap();
    group.bench_function("compiled_v3_dispatch", |b| {
        b.iter(|| compiled_v3.extract(black_box(&r1), black_box(&r2)));
    });

    // Compiled v3 via direct variant call (no enum branch)
    let simple_v3 = match &compiled_v3 {
        CompiledGeom::Simple(ext) => ext,
        _ => panic!("expected Simple variant for v3"),
    };
    group.bench_function("compiled_v3_direct", |b| {
        b.iter(|| simple_v3.extract(black_box(&r1), black_box(&r2)));
    });

    // Hardcoded v2
    group.bench_function("hardcoded_v2", |b| {
        b.iter(|| hardcoded_chromium_v2_extract(black_box(&r1), black_box(&r2)));
    });

    // Compiled v2 via enum dispatch
    let geom_v2 = parse_geometry("1{b[16]u[10]x:}2{r:}").unwrap();
    let compiled_v2 = CompiledGeom::from_fragment_geom(&geom_v2).unwrap();
    group.bench_function("compiled_v2_dispatch", |b| {
        b.iter(|| compiled_v2.extract(black_box(&r1), black_box(&r2)));
    });

    // Compiled v2 via direct variant call
    let simple_v2 = match &compiled_v2 {
        CompiledGeom::Simple(ext) => ext,
        _ => panic!("expected Simple variant for v2"),
    };
    group.bench_function("compiled_v2_direct", |b| {
        b.iter(|| simple_v2.extract(black_box(&r1), black_box(&r2)));
    });

    // Compiled v3 via extract_into (reusable buffer)
    let mut buf = ExtractionBuf::new(compiled_v3.meta());
    group.bench_function("compiled_v3_into", |b| {
        b.iter(|| simple_v3.extract_into(black_box(&r1), black_box(&r2), &mut buf));
    });

    // Stateful Extractor (owns workspace, no return value construction)
    let mut extractor = Extractor::new(&compiled_v3);
    group.bench_function("extractor_v3_full", |b| {
        b.iter(|| {
            extractor.extract(black_box(&r1), black_box(&r2));
            unsafe {
                black_box(extractor.barcode(0));
                black_box(extractor.umi());
                black_box(extractor.bio_read(0));
            }
        });
    });

    // Extractor extract-only (no accessor calls)
    group.bench_function("extractor_v3_extract_only", |b| {
        b.iter(|| {
            extractor.extract(black_box(&r1), black_box(&r2));
            black_box(&extractor);
        });
    });

    group.finish();
}

fn bench_hardcoded_vs_compiled_batch(c: &mut Criterion) {
    let mut group = c.benchmark_group("hardcoded_vs_compiled_batch");
    let batch_size = 100_000usize;
    group.throughput(Throughput::Elements(batch_size as u64));

    let reads: Vec<(Vec<u8>, Vec<u8>)> = (0..batch_size)
        .map(|i| (make_read(150, i as u8), make_read(150, (i + 1) as u8)))
        .collect();

    // Hardcoded v3 batch
    group.bench_function("hardcoded_v3_100k", |b| {
        b.iter(|| {
            for (r1, r2) in &reads {
                black_box(hardcoded_chromium_v3_extract(r1, r2));
            }
        });
    });

    // Compiled v3 batch
    let geom = parse_geometry("1{b[16]u[12]x:}2{r:}").unwrap();
    let compiled = CompiledGeom::from_fragment_geom(&geom).unwrap();
    group.bench_function("compiled_v3_100k", |b| {
        b.iter(|| {
            for (r1, r2) in &reads {
                black_box(compiled.extract(r1, r2));
            }
        });
    });

    group.finish();
}

fn bench_extract_throughput(c: &mut Criterion) {
    let mut group = c.benchmark_group("throughput");

    // Simulate processing a batch of reads
    let batch_size = 10_000usize;
    group.throughput(Throughput::Elements(batch_size as u64));

    // Chromium v3
    let geom = parse_geometry("1{b[16]u[12]x:}2{r:}").unwrap();
    let compiled = CompiledGeom::from_fragment_geom(&geom).unwrap();
    let reads: Vec<(Vec<u8>, Vec<u8>)> = (0..batch_size)
        .map(|i| (make_read(150, i as u8), make_read(150, (i + 1) as u8)))
        .collect();

    group.bench_function("chromium_v3_10k", |b| {
        b.iter(|| {
            for (r1, r2) in &reads {
                black_box(compiled.extract(r1, r2));
            }
        });
    });

    // Flex v2 with anchor
    let anchor = b"TTGCTAGGACCG";
    let sample_bc = b"SAMPLEBC10";
    let geom = parse_geometry("1{b[16]u[12]x[0-3]f[TTGCTAGGACCG]s[10]x:}2{r:}").unwrap();
    let compiled_v2 = CompiledGeom::from_fragment_geom(&geom).unwrap();
    let reads_v2: Vec<(Vec<u8>, Vec<u8>)> = (0..batch_size)
        .map(|i| {
            let gap = i % 4; // vary gap length
            (make_flex_v2_r1(gap, anchor, sample_bc), make_read(150, (i + 1) as u8))
        })
        .collect();

    group.bench_function("flex_v2_anchor_10k", |b| {
        b.iter(|| {
            for (r1, r2) in &reads_v2 {
                black_box(compiled_v2.extract(r1, r2));
            }
        });
    });

    group.finish();
}

criterion_group!(
    benches,
    bench_parse,
    bench_compile,
    bench_extract_fixed,
    bench_extract_anchor,
    bench_hardcoded_vs_compiled,
    bench_hardcoded_vs_compiled_batch,
    bench_extract_throughput,
);
criterion_main!(benches);
