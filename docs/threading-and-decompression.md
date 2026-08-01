# Threading, decompression, and the tiny-dict prefilter

Design notes for the mapping pipeline's use of threads and its k-mer prefilter.
Everything here is stated with the measurement that justifies it; where a
plausible-sounding alternative was tried and lost, that is recorded too, so it
is not retried from first principles.

Reference workloads used throughout:

| name | index | k-mers | reads | gzip structure |
|---|---|---|---|---|
| **Flex** | 10x probe panel | 1.5 M | 150 bp | multi-member, ~186 KB |
| **PBMC** | gencode v49 | 96.3 M | 91 bp | multi-member, ~12 KB |
| **bulk** | gencode v49 | 96.3 M | 150 bp | single-member |

Machine: AMD EPYC 9575F, 2 sockets x 64 cores.

## 1. Decompression is usually the constraint, not mapping

paraseq holds one mutex per input file, so exactly one thread can inflate a
given file at a time — a hard ceiling of two inflate streams for a paired run,
regardless of `-t`. On the full Flex archive (2 x 149.5 GB, 2.30 B read pairs)
that capped the run at ~11 min with mapping threads ~80% idle, and doubling
`-t` from 16 to 32 changed wall clock by 1%.

`rapidgzip-core` (feature `rapidgzip`) decodes a *single* gzip member in
parallel, lifting the ceiling without changing paraseq's structure:
**11.28 min -> 1.86 min**, output bit-identical.

Non-gzip input (plain, zstd, bz2) keeps the niffler path and allocates no
decoder threads.

## 2. Thread budget: how `-t` is split

`plan_thread_budget` splits the user's budget between mapping and decoding, and
`io::decode_budget` supervises the decoders under one shared ceiling. Upstream
adapts *within* a decoder; nothing bounds the sum *across* the 2N decoders a
multi-file run opens, and that gap is what this owns.

**Allocation matters more than budget size.** At a fixed 64 threads, the split
was worth **5.75x** (11.04 min with everything on mapping and no parallel
decode, versus 1.92 min at 32 mapping / 32 decode). Doubling 64 -> 128 threads
bought only 9%.

**Start low and ratchet.** Decoders start at 2 workers/file and double on
evidence of starvation, bounded by the budget, giving threads back only on
`Converged` (where the decoder has *measured* its own knee). Justified by an
asymmetry: past its knee every workload's curve is flat, so surplus workers
cost threads but no time, while under-provisioning costs 22-48%.

Two designs were tried first and rejected:

- **Requiring N consecutive identical verdicts.** A healthy pipeline alternates
  between `DecoderBound` and `ConsumerBound` sample to sample, so a run-length
  rule either never fires or, shortened enough to fire, oscillates. Replaced by
  an evidence score.
- **Shrinking on sustained `ConsumerBound`.** Oscillates with a ~1 s period: at
  4 workers the decoder is `DecoderBound` so the share grows to 8; at 8 it gets
  ahead of the parser and reports `ConsumerBound` so the share halves; forever.
  `ConsumerBound` is not waste — it means the decoder is comfortably ahead,
  which is the goal. Growth is therefore one-way within a run.

Knees are **not portable**: 16 workers/file for Flex, 4 for bulk, ~0 for PBMC,
because a 96 M-k-mer index makes mapping so much costlier per read that decode
demand nearly vanishes. Hence evidence-driven growth rather than a fixed ratio.

Measured (`-t 32`): Flex 111.74 s, PBMC 17.48 s, bulk 87.09 s — within 1-3% of
the best hand-tuned setting on each, with no tuning knob exposed.

## 3. The tiny-dict prefilter

A blocked Bloom filter in front of the hashbrown probe, gated on two properties.
See `tiny-dict/src/prefilter.rs` and the `FilterGate` in `tiny-dict/src/lib.rs`.

**It is consulted only on genuine searches.** `try_extend_fw` (within-unitig
extension) and piscem's fast-positive skip check both resolve against the SPSS
and return before it. So the governing statistic is the *search*-miss rate,
which differs sharply from the per-k-mer rate because extensions absorb hit runs:
Flex 99.0% (vs 79.12% per k-mer), PBMC 86.6% (vs 45.50%).

**Gate A — cache residency (8 MB).** Not a memory bound. Benefit tracks index
size, not miss rate: -16.3% CPU on the 1.5 M-k-mer panel versus -5.1%/-3.8% on
the 96.3 M-k-mer transcriptome despite comparable search-miss rates.

**Gate B — observed search-miss rate.** Sample 200 K searches, enable at >= 70%.
Gate A is the real discriminator (PBMC's 86.6% would pass B); B covers the case
A cannot see — a small index queried by reads that mostly map.

### Do not "improve" the filter without measuring first

Three changes were made and **all reverted** after measurement:

1. Multiply-shift range reduction replacing power-of-two masking. Correct in
   itself (masking inflated an 8 bits/key request to 11.2 allocated, ~40%
   waste) — but the waste was silently buying FPR.
2. 512-bit cache-line blocks with 5 bits/entry, fixing FPR 7.39% -> 2.34%.
3. Parallel construction (3.118 s -> 0.126 s at 96.3 M keys), which then
   justified relaxing Gate A.

Result: **strictly worse on every workload** — flex +8.3% CPU (from -16.3%),
pbmc +16.6%, bulk +6.2%. Both mistakes were the same shape, optimising a
quantity that was not the constraint:

- **FPR is nearly irrelevant here.** At 99% search-miss a miss is rejected
  either way, and a false positive costs only the probe that would have
  happened anyway. What is paid on *every* search is probe cost, and the
  rewrite tripled it.
- **Gate A was never about construction time.** It is about cache residency.

### External crates are not faster here

Benchmarked against `fastbloom` 0.17 ("the fastest Bloom filter in Rust") on the
real access pattern — 1.5 M keys, 99% absent queries — with a fast hasher
(rapidhash) given to both:

| | memory | build | probe | FPR |
|---|---|---|---|---|
| ours | 2.1 MB | 0.003 s | **1.68 ns/op** | **1.388%** |
| fastbloom | 1.5 MB | 0.007 s | 11.46 ns/op | 2.159% |

Not a knock on the crate — it is generic (`impl Hash` dispatch, `BuildHasher`
boundary, runtime-configured sizing). Our path is a specialised leaf: keys are
already `u64`, the hash is one `mulx`, the probe is one load and a branchless
compare. Profiling puts the filter probe at 19.4% of cycles, so the ~10 ns/op
difference is worth roughly 2% of whole-run time.

Xor/binary-fuse filters are also the wrong trade: 3 probes where we need 1, and
peeling construction of ~10 s at 96 M keys, to improve an FPR that barely
matters.

**If the remaining 19.4% is worth attacking, the lever is fewer filter probes
(the skip machinery issuing fewer searches), not a faster probe.**

## 4. Environment overrides (measurement only)

| variable | effect |
|---|---|
| `PISCEM_RAPIDGZIP_THREADS` | pins decoder workers per file, bypassing the budget |
| `PISCEM_TINY_PREFILTER` | `0` forces the filter off, `1` forces it on |
| `PISCEM_TINY_PREFILTER_MAX_BYTES` | overrides Gate A |
| `PISCEM_TINY_PREFILTER_MIN_MISS` | overrides Gate B's threshold |
