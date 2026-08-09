# Plan: closed-loop thread scheduling between decode and mapping

*Untracked working document. Not part of the repository.*

## Context

`-t` has to be divided between threads that decode gzip and threads that map
reads. That division is worth a lot — **5.75x** at a fixed 64 threads on the
Flex archive — and until now it had to be chosen before the run read a byte.

Everything we tried to choose it well failed by a wide margin:

| policy | gencode v49, 1 file, `-t 32` | full Flex pair, `-t 64` |
|---|---|---|
| closed-form supply/demand model | **60% off** optimum | — |
| two-point probe @8 workers/file | 7% off | **42% off** |
| two-point probe @budget/2 | **38% off** | 5.5% off |
| flat 1/3 of budget | 13% off | 22% off |

No variant dominated, and a plain constant had the best worst case. That is the
same shape of failure the docs already record for the original rate model — a
principled measurement losing to a simple rule — and it is why the split work
was parked.

**What changed.** Both sides are now resizable while a run is in flight:

- `paraseq` `ThreadPool` (noamteyssier/paraseq#75) — mapping workers can be
  added and retired mid-run, converging in ~161 ms, floor of one worker.
- `rapidgzip-core` `DecoderPool` (COMBINE-lab/rapidgzip-rust `feat/shared-pool`)
  — one aggregate execution budget shared across decoders, a mutable ceiling,
  and admission decoupled from the runtime throttle so a low limit no longer
  latches `Sequential`.

So the split no longer has to be *predicted*. It can be *converged to*. This
plan replaces the up-front sizing probe with a closed-loop scheduler.

## Design: what `auto` becomes

### Signals

Two, sampled on a **configurable** interval. The default is 250 ms, chosen to
match upstream's worker-retirement hysteresis — sampling much faster mostly
re-observes the same state, and 250 ms is also exactly `rapidgzip`'s
`RETIRE_AFTER`, which makes instantaneous fields alias against their own duty
cycle. It is a default, not a constant: a run whose batches take 80 ms wants a
different cadence from one whose batches take 5 ms, and callers embedding this
in another tool will have different sampling economics again. See
[Configuration](#configuration).

**Consumer (mapping) idleness** — the fraction of mapping-thread capacity not
spent mapping, from `Σ map_nanos / (threads × sample_interval)`. This is the
measurement already built for the probe; it moves from a synthetic pre-run
pipeline into the production processors, where it is one `Instant::now()` pair
per 256 records.

**Producer (decode) saturation** — from `DecoderPoolStats` and `DecoderStats`:

- `waiting_decoders` / `pool_limited` — decoders blocked on a pool slot
- `queued_tasks` — runnable work behind the admission limit
- `Σ desired_workers` vs `worker_limit` — how much more the decoders would use

`desired_workers` is the demand signal we asked for in rapidgzip#30, and it is
what makes allocation proportional instead of doubling-and-remeasuring.

### The control law

| mapping | decode | meaning | action |
|---|---|---|---|
| idle | slot-starved | decode-bound | mapping −k, pool +k |
| busy | slots unused | mapping-bound | pool −k, mapping +k |
| busy | slot-starved | balanced | hold |
| **idle** | **not slot-starved** | **source-bound (disk/network)** | **reclaim and report** |

The fourth row is new and is the strongest single argument for the closed loop:
**neither signal can identify it alone.** A probe measuring only starvation sees
idle mapping threads and spends cores on decode, which buys nothing when the
storage is the constraint. Only the conjunction distinguishes them.

`k` is proportional to the measured shortfall (from `desired_workers`), bounded
by a step cap.

The fourth row **reclaims rather than holds**, which the prototype settled. When
supply is capped elsewhere, neither side gains throughput from a thread — but
decode slots that provably cannot be used are pure waste, and on a shared machine
holding cores for nothing is the cost that matters. Reclaiming also unwinds an
overshoot, which holding does not: an over-large step leaves those threads
allocated for the rest of the run. Near the saturation point this rests within a
slot or two of the boundary, which is the right place, since "just enough to keep
up with the source" *is* that boundary.

**Watch `desired_workers` for optimism.** If a producer reports how far demand
exceeds supply rather than how many more slots would actually *help*, the broker
overshoots and then sawtooths — measured in the prototype as a limit cycling
2 → 10 → 3 → 11 → 3 for a whole run, from a model that reported `wanted: 81`
where 1 would do. `rapidgzip` should not do this (a source-bound decoder has
nothing queued, so it reports source-bound rather than a large starvation), but
it is the first thing to check if the real integration oscillates, and a damping
factor on the step is the fallback.

### Configuration

Every tuning constant is a builder default rather than a hardcoded value. Each
one below was arrived at by fixing a specific observed failure, so the defaults
carry real information — but none of them is a law, and a caller in another tool
should be able to move them without patching the crate.

```rust
let broker = ThreadBroker::builder()
    .budget(total_threads)
    .consumer(mapper)
    .producer(decoder_pool)
    // -- all optional, defaults shown --
    .sample_interval(Duration::from_millis(250))  // upstream retirement hysteresis
    .warmup(Duration::from_millis(150))           // discard the startup ramp
    .evidence_threshold(2)                        // samples before acting
    .evidence_clamp(4)                            // bound on accumulated inertia
    .deadband(0.05)                               // idle fraction treated as balanced
    .grow_step(GrowthStep::Proportional { cap: 8 })
    .shrink_step(GrowthStep::Fixed(1))            // asymmetric on purpose
    .min_consumer_threads(1)
    .min_producer_slots_per_input(2)              // below 2 rapidgzip goes sequential
    .build()?;
```

`sample_interval` and `warmup` are the two most likely to need moving, and they
interact: the warm-up should span several intervals, and the interval should
comfortably exceed the consumer's batch service time. A validating builder
should reject `warmup < 2 * sample_interval` rather than silently produce a
broker that decides on one sample.

Two defaults are load-bearing rather than merely tuned, and the builder docs
should say so:

- **`min_producer_slots_per_input = 2`.** `rapidgzip` reads the worker limit once
  while choosing a backend and commits to `DecoderPath::Sequential` below two.
  A caller who sets this to 1 gets sequential decoding for the whole file, which
  is worse than not using a pool at all.
- **`warmup > 0`.** Every workload we measured reported ~99% starved in its first
  75 ms, because the consumer is still spawning and nothing has been processed.
  Zero warm-up makes every input look producer-bound.

For piscem-rs these are also reachable as environment overrides
(`PISCEM_BROKER_INTERVAL_MS`, `PISCEM_BROKER_WARMUP_MS`) so a measurement can be
reproduced without a rebuild, matching the existing `PISCEM_*` convention.

### Guarding against oscillation

The docs already record a shrink rule that oscillated with a ~1 s period, and
now *both* sides move, which doubles the ways to get this wrong. Four
mitigations, all of which we have already had to learn once:

1. **Evidence accumulation**, not consecutive-run agreement. A healthy pipeline
   alternates sample to sample; run-length rules either never fire or oscillate.
2. **A deadband** where neither side acts, so a balanced run stops moving.
3. **Asymmetric steps** — grow the constrained side quickly, return threads
   slowly. Under-provisioning costs 22–48%; over-provisioning costs threads.
4. **A warm-up** before the first decision. Every workload we measured reported
   ~99% starved in its first 75 ms, because `paraseq` is still spawning and
   nothing has been mapped yet. Without this every input looks decode-bound.

### Starting point

Coarse and deliberately biased toward mapping: mapping gets `N − d₀`, with
`d₀ = min(2 × files, N/8)`. Two per file is the floor that keeps `rapidgzip` off
its sequential backend; the rest is earned. Growth is cheap and safe now, so the
cost of starting low is bounded by convergence time.

### What stays

- **`forced_choice`** — non-regular input, single mapping thread. Free, logical,
  no measurement.
- **The FIFO guard.** *Simpler* now: with no probe there is no prefix read, so
  the data-loss hazard largely evaporates. The per-file downgrade in
  `open_input` stays, as does the messaging.
- **`--decoder serial`** → no pool constructed at all.
  **`--decoder parallel=N`** → pin `pool.set_worker_limit(N)` and disable the
  controller. A named number means what it says.
  **`--decoder parallel`** → pool with the controller enabled.

### What is deleted

`probe_starvation`, `ProbeKernel`, `MappingProbeKernel`, `AtacProbeKernel`,
`ProbeConfig`, `Starvation`, the deadline shim, `parallel_gain`, `supply_at`,
`best_decode_share`, `MeasuredRates`, `ParallelSupply`, and the second sizing
probe. Roughly 900 lines, all of it the part we could not defend.

The measurement machinery inside it — `BatchTimer`, the windowed sampler with
warm-up — is retained and moved, because those two pieces each fixed a real
inversion.

## Genericity: recommended, with a staged rollout

The same problem appears in salmon now, and in the Rust rewrites of cuttlefish 3
(k-mer counting) and STAR (alignment) later. Four implementations of the control
law would mean making the same four mistakes four times.

**Overhead is not a concern.** The scheduler runs on one thread at the sampling
interval — 4 Hz by default, and the cost scales with it, so even an aggressive
25 ms cadence is 40 wake-ups a second on one thread. The
per-batch measurement stays in the consumer's own hot loop — it has to, and it
is already there. The generic layer sits entirely off the hot path, so dynamic
dispatch is irrelevant.

**There is enough logic to justify it.** The controller owns: windowed sampling
with warm-up, evidence accumulation, the deadband, asymmetric steps, budget
arithmetic with per-side floors, the four-quadrant decision including
source-bound detection, convergence detection, and reporting. Every one of those
is something we got wrong at least once.

### The seam

```rust
/// The side that consumes decoded bytes: mapping, aligning, counting.
pub trait Consumer: Send + Sync {
    fn set_threads(&self, n: usize);
    fn threads(&self) -> usize;
    /// Monotonic nanoseconds spent doing real work, summed over workers.
    /// The scheduler windows this itself.
    fn busy_nanos(&self) -> u64;
}

/// The side that produces them: a decompressor pool.
pub trait Producer: Send + Sync {
    fn set_limit(&self, n: usize);
    fn limit(&self) -> usize;
    fn pressure(&self) -> ProducerPressure;
}

pub enum ProducerPressure {
    /// Work queued behind the limit; `wanted` more slots would be used.
    Starved { wanted: usize },
    /// Keeping up with the consumer.
    Satisfied,
    /// Idle, but not because of the limit — the source is the constraint.
    SourceBound,
    /// Cannot use more slots at any price (e.g. a sequential-only stream).
    Inelastic,
}
```

`Inelastic` earns its place: a file `rapidgzip` decodes sequentially is 100%
starved and cannot be helped, and without this the scheduler would spend its
whole budget chasing it.

**Staging.** Build it as a crate in the piscem-rs workspace from the start
(`crates/thread-broker`), so the seam is real rather than notional, but do not
publish it until salmon consumes it. That gets the boundary right without
designing a public API for hypothetical users. Extract to its own repository
when the second consumer lands and the API has survived one migration.

## Implementation phases

Each phase ends green — tests pass, output byte-identical — so the work can stop
at any boundary.

1. **Dependencies.** Point at the rapidgzip `feat/shared-pool` branch alongside
   the paraseq fork. Confirm byte-identical output before touching anything.
2. **`crates/thread-broker`.** Traits, builder, controller, unit tests against
   fake `Consumer`/`Producer` implementations — the control law is fully
   testable without any I/O, including at compressed time scales, since a fake
   consumer can be driven far faster than 4 Hz.
3. **Consumer instrumentation.** `busy_nanos` in the three production
   processors, via the retained `BatchTimer`.
4. **Adapters.** `paraseq::ThreadPool` → `Consumer`; `DecoderPool` → `Producer`.
5. **Wire `auto`.** Replace `choose_decoder`'s probe path; keep forcings and
   preferences. Delete the probe.
6. **Validation.** Below.
7. **Docs**, then salmon.

## Validation strategy

### Correctness gates (must hold before any performance claim)

Threads moving mid-run must not change a single output byte.

| gate | how |
|---|---|
| deterministic output | RAD byte-identical at `-t 1`, every modality, both dict backends |
| output under churn | `-t 32` output identical to `-t 1` in *record content*, and record counts exact |
| no records lost | mapped + unmapped == input records, every modality |
| adversarial churn | force the scheduler to oscillate (env override driving the split 1↔N every 100 ms); counts must stay exact |
| existing parity | `--features parity-test` suites still pass |

The adversarial churn test is the important one: it is the direct analogue of
paraseq's `thread_churn_processes_every_record_exactly_once`, at the piscem
level, and it is what catches a worker retiring mid-chunk without flushing.

### Per-modality coverage

Every modality gets correctness, convergence, and optimality. Fixture status is
recorded honestly — one is missing.

| modality | index | reads | status |
|---|---|---|---|
| bulk SE | gencode v49 (96.3 M) | `human_reads/SRR21186103_1.fastq.gz` (879 MB, MarkerWindow, 88.8% mapped) | ready |
| bulk PE | gencode v49 | `human_reads/SRR21186103_{1,2}.fastq.gz` (26.1 M pairs, MarkerWindow, 70.3% mapped) | ready |
| scRNA | Flex panel (1.5 M) | `mm2/n8_R1_0.gz` + `n8_R2_0.gz`, `-g chromium_v3` (98.1% mapped) | ready |
| Flex (flagship) | Flex panel | `bigflex/big_R{1,2}.fq.gz` (150 M pairs, 2 × 9.2 GB) | ready |
| scATAC | chr1+chr2 k=23 | `atac_pbmc_5k` L001 R1/R2/R3 (3 files, 107.8 M reads) | ready |

**Bulk PE fixture — resolved.** `SRR21186103_2.fastq.gz` fetched from ENA
(877,222,851 bytes, exact match to the reported size; the local `_1` is likewise
byte-identical to ENA's). Verified: `gzip -t` clean, zero `Z_SYNC_FLUSH` markers,
`DecoderPath::MarkerWindow` at 7343 MiB/s, read names pair (`@SRR21186103.N N/1`
against `N/2`), and a full run maps 18,386,039 / 26,135,185 pairs (70.3%) in
6.14 s of mapping at `-t 32`.

This matters beyond filling a gap: bulk PE is the one modality that measures
**10.0% consumer idle** and lands on *serial*. Every other fixture is either
strongly decode-bound (Flex, 93.8%) or moderately so (bulk SE, 50.7%). Without a
mapping-bound paired workload the scheduler could converge to "always add decode
threads" and pass everything.

It also produced a result worth recording as evidence for this whole rework: the
sizing probe measured supply *falling* from 1.384 to 1.148 GB/s when given 16
workers per file, because 32 mapping threads plus 32 decode workers against a
32-thread budget oversubscribes the machine. The probe was measuring its own
interference. A closed loop never has to construct that situation, because it
only ever moves threads *between* the two sides of a fixed budget.

**Encoder check for any new fixture.** `pigz`-written gzip was unparallelizable
before rapidgzip 0.2.1 and every measurement on such a file was void. Confirm
`DecoderPath != Sequential` before using any new read set:

```bash
head -c 50000000 reads.fq.gz | python3 -c \
  "import sys; print(sys.stdin.buffer.read().count(b'\x00\x00\xff\xff'))"
```

### Thread counts

`-t 8`, `-t 32`, `-t 64` for each modality. The split behaves differently at
each: at 8 the floors dominate, at 64 the Flex case is strongly decode-bound.

## Performance measurement protocol

Three traps have already produced wrong conclusions in this work. The protocol
exists to prevent them recurring.

1. **Measure mapping seconds, not wall clock.** Index load is 1.58 s on gencode
   and is constant across configurations; including it compressed a 60% error
   into an apparent 33% one. Parse `+ Xs mapping` from the run's own report.
2. **One binary, many configurations.** There is a **~1–3% systematic
   between-binary noise floor** from code layout alone, repeatable so repetition
   does not average it away. All pinned points must come from the *same* binary
   via env override.
3. **Reach steady state.** A short run measures the startup ramp, not
   throughput. At 64 threads a sub-second run reported the *opposite* batch-size
   trend from the true one. Every configuration must run long enough that the
   ramp is a small fraction.

### The oracle sweep

For each (modality × thread count), sweep pinned decode threads across a range
spanning the expected optimum, best of 3, and record the minimum. That is the
oracle. Then run `auto` and report **distance from oracle**, which is the number
that matters — not raw time, which varies by fixture.

```
for d in <sweep>: PISCEM_DECODE_PIN=$d  -> mapping seconds, best of 3
oracle = min over d
auto   -> mapping seconds, best of 3
report (auto - oracle) / oracle
```

**Success criterion: within 10% of oracle on every modality and thread count,
and never worse than the flat-1/3 constant.** That second clause is the honest
bar — 1/3 was 13–22% off, and a scheduler that cannot beat a constant does not
justify its complexity. If it fails, say so and ship the constant.

### Convergence

Log the split trajectory each sample. For each run report:

- **time to converge** — first sample after which the split stays within ±1 for
  10 consecutive samples
- **whether it converges at all** — a run that never settles is a failure even
  if its average throughput is fine
- **final split vs oracle split**

Convergence time is the deciding factor for short runs. paraseq converges in
~161 ms and rapidgzip retires within 250 ms, so a few control steps is ~1 s. If
convergence exceeds ~20% of typical run length, that is the argument for a warm
start — and it should be made with this measurement, not assumed.

### Overhead

Two comparisons, both same-binary:

- **scheduler off vs on, same fixed split** — isolates the control loop's own
  cost. Expected ~0 at the 250 ms default. Also measure at a deliberately
  aggressive `sample_interval` (25 ms) to confirm the cost scales as expected and
  to bound what a caller can do to themselves by tuning it down.
- **`busy_nanos` instrumentation on vs off** — isolates the per-batch
  measurement in the hot path. Expect < 1%; measured at ~2 `Instant::now()` per
  256 records in the probe.

### Regression

Against current `main` (published 0.6.4 behaviour) on every modality: wall,
mapping seconds, peak RSS, and total CPU seconds. **CPU seconds matter
independently of wall clock** — the whole point of the split is to not waste
cores, and a scheduler that matches wall time while burning more CPU has failed
on a shared node.

### Robustness

| case | expectation |
|---|---|
| `-t 1` | sequential path, no pool, unchanged |
| FIFO input (`-r <(zcat ...)`) | serial for that file, no hang, other files unaffected |
| mixed regular + FIFO | per-file downgrade, correct messaging by mode |
| `Inelastic` stream | scheduler stops growing decode, warns once, does not thrash |
| cold page cache | source-bound detected; decode threads *not* grown |
| `--decoder parallel=N` | pinned, controller disabled |

The cold-cache case is how the fourth quadrant gets validated. Drop caches (or
use a file larger than RAM) and confirm the scheduler holds rather than chasing
idle mapping threads.

## Risks

**Oscillation.** The main one, mitigated four ways above, and directly tested by
the convergence criterion — a run that never settles fails the gate.

**Approximate budget arithmetic.** The pool controls *execution slots*, not OS
threads, and coordinator/scanner threads sit outside it. `mapping + pool_limit ≤
N` is a target, not a guarantee. Document it; do not pretend otherwise.

**Short runs.** A closed loop needs time. Below ~1 s it may not converge, though
the probe it replaces cost 215–500 ms on those same runs, so this is plausibly a
wash. Measure it.

**Two unreleased dependencies.** Both paraseq and rapidgzip are on branches. A
piscem-rs release is blocked until both land. Acceptable while developing;
worth stating so it is not discovered at release time.

**Genericity too early.** Mitigated by keeping the crate in-workspace and
unpublished until salmon uses it.

## Open questions

- Should the scheduler expose its trajectory in `map_info.json` for post-hoc
  diagnosis? Cheap, and it would have shortened several of the investigations
  in this work.
- Does cuttlefish 3's k-mer counting have a resizable consumer pool, or would it
  need the same paraseq change? Affects whether `Consumer` is the right seam.
