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

The bulk set was ERR3239276, which **maps at only 11%** against gencode v49
(measured at two different offsets in the file, so not a head artifact). It is
fine for I/O and threading work, where only the byte volume matters, but it is
the wrong choice for anything that depends on lookup outcomes. Use
`SRR21186103` instead: 150 bp, **88.8%** mapping rate on the same index.

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

## 2. Choosing the decoder: threads per file, not file count

The serial (niffler/libz-rs) path opens exactly one inflate stream per input
file, so its supply is `files x per-stream rate` and does **not** respond to
`-t`. Demand is `map_threads x per-thread consumption`. The parallel decoder
therefore pays off on the *ratio*, and a file-count threshold is only right at
the thread count it was calibrated at.

Measured on the 10x Flex probe panel (multi-member archives), wall-clock
speedup of the parallel decoder over the serial one:

| files | threads/file | speedup |
|------:|-------------:|--------:|
| 2 | 4 | 1.09x |
| 2 | 8 | 1.92x |
| 2 | 16 | 3.05x |
| 4 | 4 | 1.07x |
| 4 | 8 | 1.54x |
| 4 | 16 | 2.08x |
| 8 | 4 | 0.92x |
| 8 | 8 | 1.20x |

The serial path's wall clock floors at a level set purely by file count -- ~9.8 s
at one pair, ~5.0 s at two, ~3.0 s at four, regardless of `-t` -- and the
parallel decoder wins exactly where threads push demand past that floor.

The previous rule (`files >= 8` -> serial) switched on at 8 files / 32 threads
where the parallel decoder *loses* (0.92x), and off at 8 files / 64 threads
where it wins (1.20x). Replaced by `map_threads / files >= 4`. Verified:

| pairs | -t | ratio | serial | old rule | new rule |
|---|---|---|---|---|---|
| 2 | 32 | 8 | 5.05 | 3.22 | 3.22 |
| 2 | 64 | 16 | 5.21 | 2.53 | 2.51 |
| 4 | 16 | 2 | 5.23 | 5.23 | 5.27 |
| 4 | 64 | 8 | 2.95 | 2.93 | **2.55** |

**The threshold is set low (4) because the payoff is asymmetric.** Enabling the
parallel decoder when it is not needed is cheap: on a 96.3 M-k-mer
transcriptome, where mapping is slow enough that decode never binds, it measured
-4.9% wall and +2.1% CPU. Not enabling it when it *is* needed costs up to 3x
wall. The supervisor absorbs the rest -- decoders start at one worker and grow
only on demonstrated starvation, so an unnecessary parallel path degenerates to
roughly serial cost rather than to its ceiling.

**Mapping cost is the remaining term, and it is not modelled.** Per-thread
consumption ranged 0.064 GB/s (96.3 M-k-mer transcriptome) to 0.43 GB/s (1.5 M
probe panel), so the *correct* ratio is index-dependent: about 3.5 threads/file
for cheap mapping, about 37 for expensive. The ratio rule uses the cheap-mapping
figure and accepts the small loss on expensive indices, per the asymmetry above.

### Non-regular inputs (FIFOs, process substitution)

`rapidgzip` gates its parallel path on `file_type().is_file()` and falls back to
*sequential* decoding otherwise, so it has nothing to offer on a pipe. Such
inputs take the serial path unconditionally (`calibrate::Reason::NonSeekableInput`).

This also fixes a hang. `open_gz_rapidgzip` sniffs the gzip magic by opening the
path, reading two bytes, and closing it, then calls `decoder.open(path)` which
opens it *again*. On a FIFO the sniff consumes those bytes and closes the read
end, and the re-open blocks forever waiting for a writer that has exited — so
`-r <(zcat ...)` hung indefinitely with the `rapidgzip` feature enabled, while
the serial build handled it fine.

Confirmed to be ours and not upstream: `rapidgzip-core` alone, opening the same
FIFO with no prior open of the path, decompresses it correctly and byte-for-byte
identically to the regular file. Nothing was reported upstream.

## 3. Thread budget: how `-t` is split

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

## 4. The tiny-dict prefilter

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

## 5. A negative prefilter does not pay on the sshash backend

Tried and **reverted**. The motivation looked strong: dictionary seeds are
overwhelmingly misses, because a hit is absorbed into a within-unitig extension
run while every miss costs its own seed.

| workload | mapping rate | seeds / 500K reads | seed miss rate |
|---|---|---|---|
| 10x Flex R1 | 97.2% | 46.3 M | 97.7% |
| ERR3239276 | 11.1% | 53.7 M | 96.9% |
| SRR21186103 | 88.8% | 10.0 M | 82.2% |

A blocked Bloom filter over the k-mer set, probed before the seed chain, was
built and swept over 1/2/4/8/16 bits per key. Output stayed byte-identical at
every size. It lost everywhere, badly:

| bits/key | filter size | SRR21186103 CPU | Flex R1 CPU |
|---|---|---|---|
| off | -- | **30.52 s** | **48.82 s** |
| 1 | 16 MB | -- | 77.73 s |
| 2 | 32 MB | 47.85 s | 87.59 s |
| 8 | 128 MB | 43.15 s | 78.65 s |
| 16 | 256 MB | 47.87 s | -- |

Construction was not the cost -- it parallelizes to 0.09-0.12 s. Three reasons
it fails, and they compound:

- **The bucket cache already removed the work it was meant to avoid.** 71% of
  seeds now resolve their bucket from a cached value. The filter adds a probe in
  front of a chain that is no longer expensive.
- **Hashing destroys locality.** Consecutive seeds within a read touch nearby
  unitig positions, so the offsets and string reads have real locality. A Bloom
  probe is uniformly random by construction, so it is a guaranteed miss with no
  reuse -- it *adds* an access pattern strictly worse than the one it replaces.
- **Gate A already said so.** The tiny-dict filter caps at 8 MB on cache
  residency grounds. Over 96.3 M k-mers even 1 bit/key is 16 MB, so that gate
  would have rejected every configuration here. The prefilter works on a 1.5 M
  k-mer probe panel (2 MB, L2-resident) and does not generalize upward.

The lesson generalizes past Bloom filters: on this path, adding *any*
hash-addressed side structure trades a local access for a random one, and the
existing chain is local enough that the trade loses.

## 6. Sub-1% A/B differences between binaries are not trustworthy

Two separately-built binaries differ in code and data layout, and that shows up
as a **systematic** offset -- stable across runs, so extra repetitions do not
average it away and it reads exactly like a real effect.

Measured directly. Two builds differing only in the sshash SPSS decode, a
function the tiny backend never calls, were compared on the k=23 probe panel
(150M read pairs, `--dict auto` -> tiny, so the differing code cannot execute):

| | mean CPU | pairs lost |
|---|---|---|
| build A | 263.62 s | -- |
| build B | 265.54 s (+0.73%) | 4 of 5 |

An earlier round of the same comparison gave +0.77%, 3 of 4. So the floor for
between-binary comparison on this workload is around **1%**, and it is
directional and repeatable.

Within a single binary, repetition is fine -- the spread is ~0.5% -- so this is
a layout artifact, not run-to-run noise.

Consequences for reading anything else here:

- Effects **above ~3%** are safe to attribute (bucket memoization at -11.6%,
  the k=31 spill fold at -3.2%).
- Effects **near 1%** are not attributable to the code under test without
  building each variant several times with perturbed layout. A tiny-dict decode
  change measured -1.2% on gencode and +1.4% on Flex; both are inside this floor
  and neither number survives. It was reverted (`a0640b8`) on the grounds that
  no benefit was demonstrated, not because a regression was proven.
- The larger input helps precision but **not** this bias: it shrinks run-to-run
  spread while leaving the layout offset untouched.

Use `-t 1` byte-identical RAD output to confirm correctness, and reserve
performance claims for effects several times the floor.

## 7. Environment overrides (measurement only)

| variable | effect |
|---|---|
| `PISCEM_RAPIDGZIP_THREADS` | pins decoder workers per file, bypassing the budget |
| `PISCEM_TINY_PREFILTER` | `0` forces the filter off, `1` forces it on |
| `PISCEM_TINY_PREFILTER_MAX_BYTES` | overrides Gate A |
| `PISCEM_TINY_PREFILTER_MIN_MISS` | overrides Gate B's threshold |
