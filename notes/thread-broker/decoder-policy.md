# When should a run use rapidgzip at all?

Note for the thread-broker owner. Measurements plus a proposal. Nothing here is
implemented.

**Summary.** On every workload measured, choosing the parallel decoder when it
cannot add decode concurrency costs **8–28%**. Choosing serial when the parallel
decoder *could* help costs up to **4.2×**. The dividing line is mechanistic and
fits 16 of 16 measured cells, but it depends on a quantity the broker only learns
after it starts. §4 proposes making the decision reversible, which the upstream
code turns out to be most of the way toward supporting.

---

## 1. The measurements

All runs on the current `dev` binary, `--features rapidgzip-cpu-accounting`.
`serial` is `--decoder serial` (libz-rs inline, all threads map *and* decode);
`auto` is the broker.

### Flex, canonical fixture, 1 pair (100 M reads, 8.8 GB/mate)

| `-t` | serial | auto | ratio | decode slots | cost share |
|---:|---:|---:|---:|---:|---:|
| 8 | 22.29 s | 27.58 s | **1.24** | 2 | 29% |
| 32 | 22.89 s | 7.31 s | **0.32** | 9 | 28% |
| 64 | 22.51 s | 5.37 s | **0.24** | 16 | 26% |

Serial is **flat at ~22.5 s at every thread count** — one pair gives two inflate
streams, so 8, 32 and 64 threads all finish together. This is the most
decode-bound workload in the repo, and rapidgzip is worth 4.2× at `-t 64`. It
still *loses by 24%* at `-t 8`.

### Flex, `mm2` split at constant total work (4.57 GB across 1/2/4/8 pairs)

| pairs | `-t 8` | `-t 32` | `-t 64` |
|---:|---|---|---|
| 1 | 10.09 / 11.11 · **1.10** | 9.66 / 3.37 · **0.35** | 9.25 / 2.68 · **0.29** |
| 2 | 10.00 / 12.26 · **1.23** | 4.85 / 3.73 · **0.77** | 4.88 / 2.19 · **0.45** |
| 4 | 9.86 / 11.78 · **1.20** | 2.92 / 3.50 · **1.20** | 2.73 / 2.43 · **0.89** |
| 8 | 9.96 / 12.00 · **1.21** | 2.88 / 3.68 · **1.28** | 1.88 / 2.21 · **1.17** |

Serial time is almost exactly `max(D/F, M(T))` with single-stream decode
`D ≈ 9.3 s`: `9.25, 4.88, 2.73` for `F = 1, 2, 4`, then flooring on mapping.

### Overhead when the broker lands on one decode slot anyway

bulk-SE, one file, `-t 8`, mapping-bound, 3 reps:

| | mapping | vs serial | CPU |
|---|---:|---:|---:|
| serial (8 mappers) | 60.84 s | — | 474.8 s |
| rapidgzip pinned to 1 slot (7 mappers) | 64.09 s | +5.3% | 486.7 s |
| `auto` (settles at 1 slot) | 65.73 s | +8.0% | 495.8 s |

So rapidgzip's own overhead at one slot is ~5%. The rest is the partition.

---

## 2. Why serial wins when it isn't concurrency-capped

Not because rapidgzip is slow. Because **serial is work-conserving and a thread
partition is not.**

A serial mapping thread inflates under its own file's lock and then maps the
records it just produced. It is never idle: whenever there is less decode work it
does more mapping, and vice versa, rebalancing at batch granularity with no
coordination. A dedicated decode slot can only decode; when decode is momentarily
cheap it idles while mapping is short-handed, and the reverse.

So serial gives you `F` concurrent inflate streams **for free**, where `F` is the
number of decode-capable files, and costs nothing when decode is light.
rapidgzip gives you `d` streams but *takes* `d` mapping threads.

**It pays exactly when `d > F`.** Against every cell measured:

| workload | `-t` | files `F` | slots `d` | `d > F` | measured |
|---|---:|---:|---:|:--:|---|
| Flex 1 pair | 8 | 2 | 2 | no | serial 1.24 ✓ |
| Flex 1 pair | 32 | 2 | 9 | yes | rapidgzip 0.32 ✓ |
| Flex 1 pair | 64 | 2 | 16 | yes | rapidgzip 0.24 ✓ |
| Flex mm2 ×1,2,4,8 | 8 | 2–16 | 2 | no | serial 1.10–1.23 ✓ |
| Flex mm2 ×1,2 | 32 | 2,4 | 7,8 | yes | rapidgzip 0.35, 0.77 ✓ |
| Flex mm2 ×4,8 | 32 | 8,16 | 7 | no | serial 1.20, 1.28 ✓ |
| Flex mm2 ×1,2,4 | 64 | 2,4,8 | 11,13,16 | yes | rapidgzip 0.29, 0.45, 0.89 ✓ |
| Flex mm2 ×8 | 64 | 16 | 13 | no | serial 1.17 ✓ |
| bulk-SE | 8 | 1 | 1 | no | serial 1.08 ✓ |

16 of 16. This is a mechanism, not a fitted constant — which matters after the
`min_producer_slots` episode.

---

## 3. The stopgap for 0.9

`d` is what the broker *converges to*, so it is not available at open time. The
usable proxy is the budget: **require at least 8 threads per decode-capable
file**, i.e. serial unless `T >= 8F`.

That reproduces all 16 verdicts except bulk-SE at `-t 8`, which sits exactly on
the boundary and takes an 8% loss; making it strict (`T > 8F`) fixes that and
instead forgoes an 11% gain at Flex `mm2`×4 `-t 64`. Both boundary errors are
~10%, noise against the 4.2× upside.

It is calibrated on Flex, the most decode-bound workload available, so it errs
toward serial for anything less decode-heavy — the direction that protects the
mapping-bound baseline. Record the provenance next to the constant (Flex/human,
`T ∈ {8,32,64}`, `F ∈ {1..16}` files, 16 cells) **and the `d > F` mechanism it
approximates**, so the next person re-derives rather than re-guesses.

Practical effect: one paired sample needs `-t >= 16` before rapidgzip engages,
and `-t 8` becomes serial for everything. Given `-t 8` lost in all seven cells
across both workloads, that looks correct rather than merely cautious.

---

## 4. The better answer: make the decision reversible

A heuristic threshold is a prediction, and this project's history with predicting
this quantity is poor. The decision does not have to be predicted — it can be
*measured and undone*.

### The signal already exists

The broker converges to `d` within a few seconds and reports it. `d <= F` is
precisely the mechanistic statement that the parallel decoder is not buying
concurrency. So: open parallel, let the model settle, and if `d <= F`, transition
to serial and return every thread to mapping. No threshold, no workload
calibration, and the cost of a wrong start shrinks from *the whole run* to
*the detection window* — roughly 1–3 s, against a 24% worst-case penalty rate,
so a few percent overall.

### Feasibility: better than expected

Upstream already has the hard part. In `backend.rs`:

```rust
/// Resumable form of the authoritative sequential decoder.
///
/// Unlike the parallel reader coordinator, this value owns no thread. …
pub(crate) struct SequentialDecoder<C> { … }
    pub(crate) fn new(cursor, config, total_output, member_count, …)
```

It is constructed **from a mid-stream cursor plus accumulated output and member
counts**, and owns no thread. Mid-stream sequential resumption is therefore not a
new capability — the parallel path already falls back into it at eight call
sites. What is missing is only that the fallback is *involuntary*. (Historically
it was also irreversible, which is the bug that yanked three piscem releases;
that has since been fixed.)

What would need adding:

1. **An application-facing request**, e.g. `DecoderHandle::request_sequential()`,
   or giving `set_worker_limit(0)` this meaning.
2. **A quiescence protocol.** The parallel path decodes ahead of the consumer, so
   the transition point must be a drained one: stop scheduling new tasks, let the
   consumer drain in-flight results in order — which also resolves any
   marker-encoded blocks — then capture the cursor at the last delivered byte and
   construct `SequentialDecoder` from it. Draining rather than discarding is what
   makes it free.
3. **Releasing the pool permits** so the broker can hand those threads to mapping.

The pool must stay alive while any attached decoder is still parallel, which is
fine — the transition is per-decoder, and a run where one file is a tail and the
others are small is exactly a case worth handling per-file.

Risks worth naming: a single-member stream needs the 32 KiB DEFLATE window and
partial CRC at the transition offset (multi-member and BGZF can transition at a
member boundary, where state resets); and the drain adds latency proportional to
the lookahead depth, which should be bounded.

### An alternative that may be better still

Reversibility removes the *cost of being wrong*. It does not remove serial's
underlying advantage, which is work conservation (§2). The direct attack on that
is **work donation**: when a mapping thread is about to block on an empty input
buffer, let it execute a pending decode task instead of sleeping.

That is what serial does implicitly, and it would collapse the two pools into one
work-conserving pool — no partition, no threshold, and no reversibility needed,
because there would be no wrong choice to undo. It needs an additive upstream API
(roughly `pool.try_run_one_task()`) plus a paraseq change to call it on the
would-block path. The risk is priority inversion: a mapping thread entering a
long inflate when it could have been mapping. Bounded task sizes mitigate that.

I would treat reversibility as the tractable next step and work donation as the
thing to evaluate against it, not as competing proposals — reversibility is
implementable inside rapidgzip alone, while donation touches rapidgzip, paraseq
and piscem together.

### A cheap complement, either way

The rule in §3 is per-run, but the underlying problem is per-file: `F` balanced
files behave differently from `F` files where one is 100× the rest, since the
large one becomes a serial tail while the other threads idle. Compressed size is
known at open time, so the parallel decoder could be enabled only for files large
enough to be that tail. The per-file downgrade path already exists in
`open_input_pooled` (used today for non-gzip and non-seekable inputs), so this is
a small extension rather than new machinery.

---

## 5. Suggested order

1. Ship §3's threshold for 0.9, with provenance and the mechanism recorded.
2. Prototype §4's `request_sequential()` upstream; gate it on the broker's
   `d <= F` condition and re-run this grid — success is every cell at ≤1.05.
3. Evaluate work donation against it before committing to either as permanent.
4. Consider the per-file size rule once one of the above lands.

## 6. Reproducing

Fixtures: `/scratch1/rob/flex_bench/mm2/n{1,2,4,8}_R{1,2}_*.gz` (constant 4.57 GB
total, Flex library, 86.4% mapping), `/scratch3/rob/tmp/tb-flex-100m-R{1,2}.fq.gz`,
`/scratch3/rob/tmp/tb-bulk-dense-180m.fq.gz`. Indices
`/scratch1/rob/flex_bench/human_index/index` and
`/scratch1/rob/rshash_testing/piscem_wrapper_test/gencode_v49`. Raw results under
`/scratch3/rob/tmp/q{1,2,3,4}-*` and `/scratch3/rob/tmp/qf-*`.

Note the machine showed run-to-run spread of up to ~40% on short fixtures during
this session. Every cell above is ≥1.9 s and the effects are 10%+, but a formal
gate should use paired counterbalanced controls rather than the sequential runs
used here.
