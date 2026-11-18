# Experimental SSHash minimizer fingerprints

The Cargo feature `experimental-bucket-fingerprints` eagerly enables SSHash's
derived minimizer filter when `ReferenceIndex<Dictionary, ContigTable>` loads.
It is a research switch, off by default, with no runtime flag or automatic
amortization policy. The loader logs the allocated payload bytes. Index files,
skipping algorithms, mapping options, and RAD semantics are unchanged. Independent
seed batching is not integrated into Piscem by this patch.

## Dependency and build requirements

This feature requires the unpublished `experimental-query` API from
[SSHash PR #10](https://github.com/COMBINE-lab/sshash-rs/pull/10), branch
`perf/sbwt-seed-campaign`. This is a **dependent draft change** until
that API is released and Piscem's dependency can be updated. In particular,
`--all-features` against the current registry dependency does not work. Ordinary
default-feature registry builds continue to resolve and build.

Clone/check out the matching SSHash branch, then from this repository run:

```bash
cargo build --release \
  --features experimental-bucket-fingerprints,sshash-lib/experimental-query \
  --config 'patch.crates-io.sshash-lib.path="/absolute/path/sshash-rs/crates/sshash-lib"'
```

Both feature arguments are required. We deliberately do not forward the dependency
feature from `Cargo.toml` yet: Cargo resolves inactive feature edges too, so naming
an unpublished dependency feature there would break ordinary registry builds.
No local path patch is committed. Use a separate checkout/target directory for
research builds and keep their locally resolved lockfiles out of the PR. After
publication, the intended cleanup is to update the dependency and forward the
upstream feature normally, then exercise it in CI.

## Runtime behavior and costs

```bash
# The index and reads must be the same for original and filtered runs.
piscem-rs map-bulk --dict sshash --decoder serial -t 1 \
  -i /data/index -1 /data/reads_1.fq -2 /data/reads_2.fq -o /tmp/filter-run
```

Specify `--dict sshash` for these measurements. With automatic dictionary choice,
a loaded SSHash dictionary may subsequently be replaced with Tiny, wasting the
filter's initialization work. The hook applies whenever this SSHash index loader
is used; it does not add filtering to Tiny.

The filter is a boxed byte array indexed by the existing minimizer MPHF. Light
buckets store seven-bit tags; heavy/empty buckets have a bypass marker. A mismatch
avoids bounds, anchor, and text accesses. A matching fingerprint still receives
full verification. Heavy buckets bypass rejection to preserve conservative
`minimizer_found` hints, not just membership. All original lookup fields are
preserved; the table is derived and never serialized.

Construction visits buckets sequentially and reads the first light-bucket anchor
and sequence. The measured Gencode v44 dictionary allocated 13,228,125 bytes
(12.6 MiB) and took roughly 0.21 seconds to initialize the filter. A separate human
dictionary required 177,989,834 bytes (169.7 MiB) and about 11.5 seconds. Short jobs
can lose more time in initialization than they save in mapping.

## Recorded campaign

Baseline Piscem: `9d3bd92b3611590cf3c6514a604e290008881812` (0.12.3).
Baseline SSHash: `80131cd3806e79d4d74fa7cbceff9e5b8da16a8d` (0.8.1).
The existing Gencode v44 cuttlefish input was rebuilt with the current unmodified
baseline to obtain a compatible PCTAB04 contig table and single-partition MPHF;
older PCTAB03 files cannot be loaded by this Piscem version.

The experiment used an EPYC 9555, Rust 1.97.1, thin LTO, and
`-C target-cpu=native -C llvm-args=-x86-branches-within-32B-boundaries` for both
variants. Each dataset contained 100k paired reads (simulated, real RNA, or random),
uncompressed FASTQ, one mapping worker, serial decoding, explicit SSHash, and the
default permissive skipping policy. Nine repetitions alternated variant order.

Before starting, three 10-second quiet samples required host busy <4% and benchmark
CPU/SMT sibling busy <5%. Queries were pinned to CPU 24. All 54 timed runs passed
the subsequent acceptance thresholds: <2 busy CPU-equivalents outside CPU 24 and
<=1% activity on sibling CPU 88 during overlapping telemetry samples.

| Dataset | Original mapping s | Filter mapping s | Original wall s | Filter wall s | Mapped pairs |
|---|---:|---:|---:|---:|---:|
| Simulated | 0.49297 | 0.48335 | 0.54091 | 0.74791 | 96,568 |
| Real RNA | 0.49679 | 0.45225 | 0.54464 | 0.71266 | 83,996 |
| Random | 1.32654 | 0.61736 | 1.37444 | 0.88266 | 0 |

Values are medians. Mapping time excludes filter initialization; wall time includes
it. Real mapping improved about 9%, while the whole short job was about 31% slower.
Random mapping improved about 2.15x, or 1.56x including initialization. These data
do not establish a benefit for arbitrary read sets, thread counts, or job lengths.

Every RAD file was compared byte-for-byte against the first run of its dataset,
across both variants and all repetitions. Semantic metadata (counts, signatures,
mode, k/m, reference count, thread count, and skipping policy) was also equal.
Timing metadata is deliberately excluded from equality. The retained
[per-run records](experimental-bucket-fingerprints/results.json) include timings,
counts, byte sizes, SHA-256 hashes, and acceptance flags. Multi-worker output byte
order was not tested by this campaign. The feature-enabled library suite also
passed 283 tests with three ignored during the campaign.

To repeat on another dataset, build original and candidate with identical flags
and compiler, retain both binaries, wait for a quiet physical core, alternate
their order, and use a distinct output prefix per run. Compare complete `.rad`
bytes with `cmp`, semantic `.map_info.json` fields separately, and report both
`mapping_seconds` and measured process wall time. Do not compare runtime fields
for equality or claim zero-read mapping from file size: RAD headers still contain
reference information. Collect fresh host/sibling telemetry for acceptance.

## Configuration and follow-up

The upstream library already provides runtime enable/disable methods. For Piscem,
a runtime off/on option and explicit memory budget would be more useful than a
compile-time switch once the dependency is published. Keep the default off until
startup amortization and workload effects are understood; no automatic policy is
implemented here.

The seven-bit layout is not proven optimal. Wider tags only reduce collisions
for absent selected minima landing in light buckets; they cannot reject an absent
k-mer whose minimizer is genuinely indexed or a heavy-bucket bypass. Upstream
documents packed nibble, 255-tag byte, and wider-entry experiments, but none are
implemented or measured in this change. Useful next measurements include filter
outcome counters in a separate diagnostic build, layout/memory tradeoffs, bounded
parallel initialization, and batching independent seeds while preserving each
read's streaming state. See the companion SSHash documentation for the direct
SBWT comparison and why substring-based skipping requires a stronger absence
certificate than this minimizer filter provides. The full
[upstream design and experiment notes](https://github.com/COMBINE-lab/sshash-rs/blob/perf/sbwt-seed-campaign/docs/experimental-query.md)
include API examples, tuning alternatives, exactness tests, and standalone harnesses.
