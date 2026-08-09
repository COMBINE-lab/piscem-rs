# Thread-broker design and implementation audit

Date: 2026-08-06
Scope: commit `ac451f7` on `feat/thread-broker`, including
`crates/thread-broker`, the piscem adapters, the three mapping CLIs, the two
resizable upstream pools, tests, and the accompanying design/validation docs.

## Executive assessment

The core idea is strong enough to continue: consumer starvation alone cannot
identify the balance point, measuring non-blocked work is the right direction,
and a model-based jump guarded by observed throughput is much better founded
than the controller it replaces. The design document is also unusually useful:
it records failed approaches, exposes assumptions, and names falsification
criteria.

The implementation is nevertheless **not release-ready**. There are confirmed
budget-plumbing defects in production paths, the claimed budget invariant is
only an invariant of requested targets rather than actual execution, actuation
failures can be silently ignored, and the producer measurement is an aliased
sample rather than the cumulative busy-time signal the model really needs. The
permanent tabu behavior also does not have the bounded-neighbor fallback claimed
by the design.

The recommended disposition is:

1. Keep the closed-form estimator and the generic crate.
2. Block release on gates A-C below: authoritative budget plumbing, actual
   budget enforcement, and explicit actuation failure handling.
3. Treat the measured 3-5% oracle gap as preliminary until producer busy time,
   allocation-dependent scaling, and the full modality/thread-count matrix have
   passed gates D-H.

## What is already good

- The diagnosis of the one-sided starvation controller is persuasive and is
  pinned by `does_not_walk_away_from_the_optimum`.
- `Work` clearly distinguishes non-blocked work from wall time and keeps the two
  stages' progress units independent.
- `BusyMeter` publishes inside a batch, avoiding the previously observed zero
  deltas on long batches.
- Warm-up, smoothing, post-move blackout, like-for-like ratification windows,
  and shrink-before-grow *request ordering* each answer a real failure mode.
- The cap has finite history instead of a permanent running minimum.
- The final `Model` is exposed, making a bad decision diagnosable rather than
  opaque.
- The controller fake couples producer supply and consumer demand; it is much
  better than an independently scripted utilization fake.
- The design already acknowledges the major open questions: steady-state flow,
  cost-share inconsistency, absolute thresholds, source-bound inference,
  modality coverage, and overhead.

Those strengths should survive the corrective work below.

## Findings

### TB-01 — Critical: there is no single authoritative effective budget

Several production paths do not divide the same `N` that the user or resource
environment selected.

1. `plan_thread_budget` caps `total` with `available_parallelism` at
   `src/io/fastx.rs:288-325`, but `ThreadBudget` does not retain that effective
   total. Bulk and scRNA subsequently size both pool maxima and the broker with
   the original `num_threads` (`src/cli/map_bulk.rs:420-467` and
   `src/cli/map_scrna.rs:518-566`). The broker can therefore expand immediately
   past the cgroup/cpuset-aware cap the planner just computed.
2. When decoder selection is serial, all three CLIs set `decode_budget = 0` but
   do not return the reserved opening slots to `map_threads`. The mapping pool
   permanently starts below the effective budget. The same under-allocation
   occurs in builds without the optional `rapidgzip` feature: a decode share is
   planned, no decode pool or broker exists, and the mapping pool remains at the
   reduced target.
3. scATAC rebinds `num_threads` to `plan.map_threads` at
   `src/cli/map_scatac.rs:436`, then uses that smaller number as both the mapping
   maximum and broker budget at lines 470-479. In adaptive mode, part of `-t` is
   permanently lost and the log reports the smaller number as though it were the
   user's budget.
4. `-t 1 --decoder parallel` reaches `clamp(1, num_threads - 1)` with an invalid
   range, and `parallel=N` also subtracts one from a one-thread budget. The CLI
   promises `-t 1` behavior, so this must not panic.
5. `PISCEM_RAPIDGZIP_THREADS` is documented in `plan_thread_budget` as an exact
   override, but it neither obeys the usable-core cap nor disables the broker.
   The broker can overwrite it, and an oversized value can make
   `map_threads + decode_budget` exceed the nominal budget.

This is the most urgent finding because every performance claim depends on all
compared runs using the same resource envelope.

**Improvement.** Introduce one immutable execution plan, for example:

```text
ExecutionPlan {
    requested_budget,
    effective_budget,
    mapping_target,
    decode_target,
    decode_mode: Serial | Adaptive | Pinned,
}
```

Construct it once after feature availability, input forcing, explicit user
preference, and usable-core detection are known. Every pool maximum, initial
target, broker budget, log message, and decoder headroom must derive from it.
Do not shadow the budget in a CLI. Decide whether a pin outside the available
budget is rejected or clamped with a warning; never silently reinterpret it.

**Acceptance gate A.** Table-driven tests over all three CLIs (or a shared pure
planner they all call), feature on/off, budgets `{1,2,8,32}`, usable-core caps
`{1,4,16}`, input outcomes `{serial, adaptive}`, preferences
`{auto, serial, parallel, parallel=N}`, and the environment override must prove:

- no panic and no invalid clamp;
- serial/no-feature: `mapping_target == effective_budget`, `decode_target == 0`;
- adaptive/pinned with `effective_budget >= 2`: both floors hold and
  `mapping_target + decode_target == effective_budget`;
- every pool maximum and `ThreadBroker::budget` equals `effective_budget`;
- the three CLIs produce the same plan for equivalent inputs;
- `-t 1` always selects a defined sequential plan;
- invalid pins are rejected or visibly clamped according to one documented
  rule.

### TB-02 — Critical: requested-target arithmetic is not actual budget enforcement

`ThreadBroker::apply` requests a shrink and immediately requests growth on the
other side (`crates/thread-broker/src/controller.rs:438-459`). Both upstream
contracts explicitly allow lag:

- paraseq workers retire only after finishing a batch;
- rapidgzip tasks admitted before a lowered limit keep their slots until the
  CPU-intensive task finishes.

Consequently, setter call ordering does **not** prevent transient
oversubscription. A 16,384-record mapping batch can make that transient material.
The test `never_oversubscribes_the_budget` cannot detect it because both fake
setters apply synchronously and it asserts only the broker's requested tuple.

There is a second upstream constraint. A paraseq `Collection` divides a pool
among concurrently active reader groups. Each share has a floor of one live
worker, so the true consumer floor is the number of active groups, not
`BrokerConfig::min_consumer_threads == 1`. If the broker requests fewer mapping
threads than active groups, the pool cannot comply. Moreover, integer division
means the aggregate live target is generally not exactly the requested target.
The broker reads aggregate live workers for diagnostics but does not use that
fact to constrain or acknowledge a move.

Finally, decoder coordinator/scanner threads and the sampler are outside the
controlled slot sum. That can be a documented distinction, but then `-t` is an
execution-slot budget rather than a literal thread/CPU budget and all performance
comparisons must account for the extra CPU.

**Improvement.** Make resizing a two-phase state transition:

1. request the shrinking side;
2. observe an acknowledgement or actual active count at/below its new target;
3. only then grow the other side;
4. time out as an explicit actuation error rather than silently proceeding.

Expose the consumer's true dynamic floor (at minimum, active reader groups) and
the producer's actual occupied slots through the broker contract. Alternatively,
fix the upstream pool to enforce a genuinely global aggregate target, but it
must still preserve at least one worker for every input group that is expected to
make progress.

**Acceptance gate B.** Add event-level instrumentation and an adversarial real-
pool test that drives `1 <-> N-1` every 100 ms with `{1,2,8,24}` logical input
groups. Across every resize event:

- `actual_mapping_live + actual_decode_busy <= effective_budget` at all times;
- neither side falls below its declared dynamic floor;
- every requested target is either acknowledged within a measured timeout or
  returns a surfaced error;
- every record is processed exactly once and output content matches a fixed
  run;
- steady-state process CPU concurrency is no more than
  `N + max(1, 0.03*N)` unless a different auxiliary-thread allowance has been
  explicitly documented and measured.

If strict instantaneous enforcement is intentionally rejected, replace the
first criterion with a quantified oversubscription allowance (peak, p99, and
slot-seconds above budget) and stop claiming that shrink-before-grow prevents
it.

### TB-03 — High: actuation and sampler failures are silently converted into success

The generic traits return no result from `set_threads` or `set_limit`. The real
decode adapter logs and ignores a refused pool limit
(`src/io/broker.rs:139-149`), while the broker updates its internal `split` as if
the request succeeded. This is the exact divergence the adapter documentation
warns about.

`ThreadBroker::start` discards sampler thread-spawn failure with `.ok()`
(`controller.rs:151-158`). `RunningBroker::finish` converts a sampler panic into
`BrokerReport::default()` (`controller.rs:75-83`). A user can therefore receive a
zeroed report rather than an error, after the opening split has already been
applied.

**Improvement.** Make both resize methods fallible and return the accepted
target or an acknowledgement token. Make `start` and `finish` return
`Result<_, BrokerError>`. On a partial two-sided move, either complete a safe
rollback or enter an explicit degraded fixed-split state. Include requested and
observed targets in the error/report.

**Acceptance gate C.** Failure-injection fakes must refuse each half of growth,
shrink, and rollback in turn; a panic-injected sampler must also be tested. Every
case must satisfy all three properties:

- no silent default report or merely logged warning;
- broker state and observed pool state cannot diverge without an error;
- the controlled-slot budget remains valid or the job aborts before processing
  continues.

Also add a test where `min_consumer_threads > 1` and the requested initial
producer split is `budget - 1`; `start` currently clamps producer slots without
respecting the consumer floor and can request a sum above the budget.

### TB-04 — High: the permanent tabu set can make the broker permanently deaf

After a regression, the exact producer target is inserted into a permanent
`HashSet<usize>`. In `Survey`, seeing that target again immediately enters
`Steady`; in `Steady`, the same target is filtered out before it can trigger a
resurvey (`controller.rs:339-341` and `381-386`).

The design says neighboring targets remain available and bounds the damage to a
few threads. The implementation does not select a neighbor: the solver returns
one rounded target, and if that target is tabu no alternative is tried. A noisy
or transient ratification can therefore prevent the broker from ever reaching a
later regime whose true optimum is exactly that value.

**Improvement.** Scope rejection to a workload/model epoch and retain the reason,
baseline, achieved rate, and confidence. Within one epoch, choose the nearest
non-tabu candidate or the best previously measured split instead of sleeping at
an arbitrary `from`. On a statistically established regime change, expire old
rejections or revalidate them after a cooldown. This preserves termination
without permanent deafness.

**Acceptance gate D.** Add deterministic scenarios in which:

- target `T` is rejected because of a temporary post-move slowdown, then becomes
  the true optimum after a regime change;
- the biased model repeatedly proposes `T` while `T +/- 1` is safe;
- two regimes legitimately revisit the same optimum.

The broker must never retry the same failed target within one unchanged epoch,
must reconsider it after a proven epoch change, must finish within 10% of the
oracle, and must make at most five moves per epoch.

### TB-05 — High: producer busy time is not measured with known accuracy

`DecodeProducer::work` integrates an instantaneous `busy_workers` snapshot with
a left Riemann sum (`src/io/broker.rs:85-121`). At the 100 ms sampling interval,
short tasks and periodic decoder behavior can be missed or overcounted. The
controller solves from only the most recent three accepted windows, so the
design's claim that sampling error averages over the whole job does not protect
an individual decision. Choosing 100 ms rather than 250 ms also does not remove
aliasing; their phases repeat.

The signal excludes coordinator/scanner work and copies performed outside pool
slots. The real adapter test proves only that the counter is nonzero, monotone,
and below `limit * wall`; it does not compare it with ground truth. The reported
37% versus 26% cost shares are therefore credible evidence of either measurement
bias or allocation-dependent service cost, not a small calibration detail.

Consumer instrumentation also measures only the batch callback. Per-thread state
creation occurs before the timer, and `on_batch_complete`/`on_thread_complete`
RAD flushing occurs after it (`src/mapping/processors.rs:326-402` and analogous
paths). The design's claim that state construction contaminates consumer busy
time is therefore not true for the current timer placement. Whether excluded
output locking is correct depends on whether the model intends to price all
consumer occupancy or only its scalable portion; that needs an explicit rule.

**Improvement.** Prefer a cumulative upstream counter updated on pool permit
acquire/release, ideally split into scalable decode work, auxiliary CPU work,
and blocked time. Use it directly instead of reconstructing occupancy by polling.
Add parallel component timing on the consumer side long enough to quantify state
construction and output flush cost before deciding whether they belong in the
model.

**Acceptance gate E.** Against an event- or CPU-time ground truth on dense,
marker-window, sequential, bursty, mixed compressed/plain, and source-bound
inputs:

- whole-run producer busy-time relative error <= 3%;
- p95 error of the exact three-window aggregate used for a decision <= 10%;
- changing sample intervals across `{25,50,73,100,137,250}` ms and randomized
  phase changes the solved cost share by <= 2 percentage points;
- accounted producer components explain at least 95% of measured producer CPU;
- excluded consumer work is either < 2% of consumer CPU or is explicitly modeled.

If these gates fail, the upstream cumulative signal is a release prerequisite,
not an optional refinement.

### TB-06 — High: budget-invariant cost share is not a valid general expectation

The closed form assumes each stage's per-item cost is independent of the number
of threads assigned. That is plausible for an ideal linear stage, but mapping a
large shared index can become memory-bandwidth/cache-contention limited, output
has a serialized component, and decoding can have gradual diminishing returns
well before a hard source cap. In those cases `s_c` and `s_p` are functions of
allocation.

Therefore, the design's statement that the same workload's cost share *should*
agree at `-t 32` and `-t 64` is too strong. The 37%/26% difference could be an
instrumentation bug, but it could also be a real change in thread-time per record
under contention. The current controller implements a hard producer cap but no
observed scalability curve for either stage, despite citing that production
guard in the design.

**Improvement.** Measure stage service cost and end-to-end throughput at pinned
splits, not just across total budgets. If cost changes materially with
allocation, retain `(allocation, cost/rate)` observations and solve against
monotone empirical response curves. A simpler fallback is a bounded local check
around the model answer rather than trusting one allocation-independent jump.

**Long-run evidence (2026-08-07).** A randomized three-repetition scATAC sweep
at `-t 8` and 2,000,000 records found producer 5 to be the stable wall oracle
(13.662 s median). Adaptive always followed 4→6→7→6 and finished at 16.873 s,
23.5% slower. The current geometric correction cannot discover producer 5
because it never measures the interior of a retained two-slot jump. This is now
a demonstrated release failure, not a hypothetical scaling risk. The bounded
local refinement is now implemented: four long real runs all followed
4→6→7→5, converged in three moves with no resurvey, and had a paired median
mapping ratio of 1.0123 (four-pair upper 1.0501) against pin 5. A subsequent
full-calibration-freeze gate exposed a second-order experimental bias: probing
5 directly from failed endpoint 7 cold-started two consumer threads but compared
them with a warm rate retained at 6. The controller now restores and blacks out
at the retained baseline before testing its interior neighbor. The model-only
freeze policy is still unsafe for this workload: four balanced long runs moved
to producer 1, were 42.2% slower, and used 3.21x the oracle's aggregate CPU
despite consuming only 0.401 ms median administrative CPU.

`FreezeAfterFullCalibration` now closes the safe low-overhead path. Its interior
candidate uses independent, bounded confirmation horizons to distinguish
post-resize ramp from settled response without adding steady-state samples. In
the final four-pair long gate, every run selected producer 5 in four moves,
converged in 3.132--4.125 s, stopped monitoring, and matched canonical output.
The paired mapping median/upper ratios were 1.0162/1.0369 versus pin 5, and the
worst controller+sampler lifetime CPU was 2.582 ms. This satisfies the
scATAC-specific freeze acceptance gate; the model-only freeze remains available
only as the explicitly cheaper, lower-fidelity policy.

Stable scATAC now defaults responsive monitoring to 5 s while preserving a
same-binary 25 ms override. An eight-pair balanced direct A/B retained producer
5 and canonical output in every run. Sparse/normal upper ratios were 1.0211
mapping, 1.0209 process wall, and 1.0197 aggregate CPU, while median
administrative CPU fell from 4.902 to 2.623 ms and monitoring observations from
431 to 10. This passes the policy-specific 5% non-regression and direct
0.001-core gates.

**Acceptance gate F.** For each validation modality and each `N`, pin at least
five splits spanning both sides of the optimum and record producer/consumer busy
time per logical record, throughput, CPU, and buffer occupancy. After accounting
for measurement error:

- if either per-item cost varies by > 10% across the useful range, the constant-
  cost model must not ship without an allocation-aware correction;
- the corrected model's predicted throughput at each held-out split must be
  within 10% of measured throughput;
- `auto` must meet gate H even when stage scaling is deliberately sublinear in
  simulation;
- a retained nonlinear jump spanning more than one slot must locally refine
  every unmeasured interior candidate needed to identify a discrete peak. On the
  long scATAC fixture, `auto` must measure producer 5, reach a final split within
  5% of its median wall time, keep the one-sided 95% upper ratio <=1.10, and
  converge within five moves/five seconds without resurvey;
- any freeze policy offered for scATAC must complete that nonlinear refinement
  before stopping, use <=5 ms total administrative CPU, and emit no recurring
  samples after convergence.

Do not use cross-budget cost-share invariance as the instrumentation gate. Use
repeatability at the same split and agreement with cumulative ground truth.

### TB-07 — Medium-high: flow transients, fixed blackout, and the cap can reinforce one another

The busy-ratio substitution is exact only when equivalent logical work crosses
both stages over the measurement span. Multiple buffered readers, file
boundaries, and long runs of different compressibility violate that premise.
The current three-window smoother has no queue-occupancy or in-flight-work guard
to determine the size of the mismatch.

Blackout is a fixed sample count rather than an observed convergence condition.
A slow 16,384-record consumer batch can outlast the default 400 ms, while a fast
resize can settle much sooner.

The slack cap is triggered when **any one** of 32 recent windows has more than
one worker of slack (`controller.rs:689-697`). Unlike pressure saturation, it
requires no persistence. An ordinary lull or buffer transient can thus create a
multi-second cap using the same aliased producer signal. History is counted in
samples rather than time, so changing `sample_interval` silently changes cap
memory from 0.8 s to 8 s over the supported interval range above.

**Improvement.** Record buffer/in-flight deltas or form cost windows between
checkpoints with comparable queue occupancy. End blackout only after both sides
acknowledge their targets and a minimum clean window has elapsed. Make cap
evidence persistent, duration-based, and conditional on comparable grants;
retain confidence and reason in `Model`.

**Acceptance gate G.** In deterministic simulations and a real multi-file
fixture with alternating compressible/incompressible regions and 0.25-1.0 s
lulls:

- clean-window cost-share bias stays <= 5 percentage points;
- a source ceiling at `K` is identified within 1 s and caps at no more than
  `K + 1`;
- a transient lull cannot keep the producer below its later required allocation
  for more than 1 s after work resumes;
- changing `sample_interval` does not change cap retention in wall-clock time by
  more than 20%;
- no decision uses pre-acknowledgement resize data.

### TB-08 — Medium: ratification has no confidence model despite permanent consequences

Ratification compares one recent block with one post-move block and applies a
fixed 5% tolerance. Post-move windows are accumulated whether or not
`usable_window` accepts them. A single I/O pause, output flush, or regime boundary
can therefore reject a target permanently. Conversely, a bad move inside the
noise band is always kept. This is defensible as a provisional guard, but not as
strong evidence for a permanent tabu entry.

**Improvement.** Require a minimum item count and usable windows, retain block
samples, and use a robust paired/block comparison. Distinguish `regressed`,
`inconclusive`, and `kept`; only a confident regression should create tabu
state. Inconclusive moves can remain in place but should not become permanent
evidence.

**Acceptance gate.** Across at least 1,000 seeded noisy traces calibrated from
real window variance, false confident rejections of equal-throughput moves must
be < 1%, and a true >=10% regression must be detected in >=95% of traces within
1 s. Isolated zero-progress and output-flush windows must not create a tabu
entry.

### TB-09 — Medium: absolute hysteresis values do not scale, and some report semantics are misleading

The design already recognizes that `deadband_threads = 2` and
`resurvey_distance = 6` are not derived. At `-t 8`, those are 25% and 75% of the
budget. Persistence reduces noise-triggered resurveys but does not make these
distances scale appropriately.

`BrokerReport` documents trajectories as decisions, but the code appends one
entry every post-warm-up sample. `converged()` reports only that
`time_to_converge` is set; it does not expose whether the controller ended in
survey, blackout, ratification, steady, or accumulating drift. Long runs also
retain unbounded per-sample trajectories.

**Improvement.** Derive movement/reopen thresholds from uncertainty in cost
share or solved target, with absolute floors only for tiny budgets. Report final
phase, current drift duration, requested versus observed split, cap reason, last
actuation error, and peak controlled-slot occupancy. Make full traces bounded or
opt-in and expose them in a machine-readable artifact.

**Acceptance gate.** On budgets `{2,4,8,16,32,64,96}`, stable noisy simulations
must have zero false resurveys in 30 s for at least 99% of seeded runs, while a
regime shift of at least 10% of budget must reopen within 2 s in at least 99%.
The report must distinguish every terminal phase and must not label an active
survey/ratification or sustained pre-threshold drift as unqualified convergence.

### TB-10 — Medium-high: the tests validate arithmetic more than real robustness

The 16 controller tests pass, but the fake has the same constant-cost equation
as the solver, immediate resize, perfect shared progress, no queue, no burstiness,
and no partial scalability. It cannot expose TB-02, TB-05, TB-06, or TB-07.

`re_solves_after_a_regime_change` mutates `Pipeline::s_c` through a raw pointer
while the broker thread can read it. The state mutex does not protect the field
itself, so this test contains a Rust data race and its result is not valid
evidence. Move mutable costs inside the mutex or use atomics.

The real adapter tests are ignored by default and use machine-specific fixture
paths. Even the fixture-free limit round-trip test is ignored. There are no
unit tests for `plan_thread_budget`, no cross-CLI execution-plan tests, no real
multi-share consumer test, and no end-to-end broker correctness/churn test from
the validation plan.

**Improvement and gate.** Build a deterministic discrete-event pipeline fake
with queue capacity, delayed/refused resize, nonlinear rates, lulls, counter
resets, and seeded noise. Make small dense/marker/sequential gzip fixtures
portable or provision them in CI. Run the feature matrix in CI. No controller
unsafe code should be needed in tests. Gates A-H should each have at least one
test that fails against the current implementation.

### TB-11 — Medium: current documentation and measurement instructions are not reproducible

`docs/threading-and-decompression.md`, CLI help, and the top of
`src/io/calibrate.rs` still describe the deleted pre-run probe and old
supervisor. The validation plan invokes `PISCEM_DECODE_PIN`, which does not exist
in the tree. `PISCEM_RAPIDGZIP_THREADS` is documented as a pin but does not
disable the broker. `parallel=N` is described as a per-file ceiling, while the
pooled implementation multiplies it into an aggregate limit and gives each
decoder headroom up to the total budget. Non-gzip/stream inputs are included in
that multiplication even though they cannot use pool slots.

This is not cosmetic: it prevents another reviewer from reproducing the oracle
sweep and makes the public override semantics ambiguous.

`RUSTDOCFLAGS='-D warnings' cargo doc -p piscem-rs --no-deps --features
rapidgzip` currently fails with five documentation errors: four deleted-probe
links in `src/io/calibrate.rs` and one public-to-private `BusyIntegrator` link.

**Improvement.** Add a single aggregate fixed-producer control intended for
measurement, clearly separate it from any per-file decoder headroom, and make it
disable the broker. Record effective plan and broker trace in `map_info.json` or
an adjacent JSON file. Rewrite user docs and CLI help after semantics are fixed.

**Acceptance gate.** A documented local harness must reproduce every oracle point from
one binary, verify requested and actual slot counts before accepting timing, and
emit raw data plus commit/dependency revisions. Rustdoc with warnings denied and
CLI help snapshot tests must pass in both feature configurations. Both upstream
branch dependencies must be pinned to reviewed commits and either merged/released
or explicitly accepted as a release blocker.

### TB-12 — Medium: broader overhead and short-run coverage remain incomplete

The producer sampler now reads each handle's lock-free executing-worker signal
at a jittered 3 ms cadence only during calibration, then 25 ms in ordinary
monitoring; stable scATAC responsive operation probes every 5 s by default. The
in-batch meter reads the clock and updates busy time every 256 records. During an
open scATAC decision, completed progress is published every 64 records into
cache-padded per-processor shards. Each shard has one writer and therefore uses
a relaxed cumulative store; only the controller sums shards, at its sampling
cadence. This removes the shared-cache-line contention found by the first formal
fine-publication run. Jobs shorter than the warm-up plus
survey/blackout/ratification sequence may still never benefit from adaptation.

**Acceptance gate.** Same-binary, randomized fixed-split comparisons with the
broker off/on must retain <=1% one-sided upper bounds for mapping wall, process
wall, and aggregate process CPU at 100 ms for `{1,2,24}` inputs;
instrumentation on/off must also be <=1%. Aggregate process CPU is only a
regression backstop, because 1% of a saturated 64-thread process can equal 0.64
core. Direct two-read lifetime CPU accounting for the controller and sampler
must additionally show <=0.001 core (1 ms CPU per wall second) for responsive
operation, while freeze-after-convergence must consume <=5 ms administrative
CPU total and no recurring CPU afterward. At 25 ms, record the scaling rather
than requiring the same process-level gate. For jobs of `{0.25,0.5,1,2,5}` s,
`auto` must be no more than 5% slower than the safe opening split, or adaptive
mode should be bypassed/warm-started below a derived duration threshold.
The scATAC fine-publication sub-gate uses a documented local, exactly
position-balanced 30-pair fixed-split 64-versus-256 comparison with the same
three <=1% bounds and canonical output. The final 2-million-record run had
exactly 15 observations in each first-position stratum and matching canonical
output/counts in all 60 runs. Fine/generic paired medians were `0.99933` mapping
wall, `0.99895` process wall, and `0.99965` aggregate CPU; design-aware
one-sided 95% upper ratios were respectively `1.00419`, `1.00434`, and
`1.00441`. This sub-gate passes. The broader duration/input-count policy matrix
and remaining short-run cells are still required. Because a pinned run has no
controller to change measurement phase, this A/B deliberately holds the fine
cadence for the full job; the adaptive path returns to 256 after calibration.
Unless an explicit validation override is present, pinned/no-broker execution
should therefore select 256 immediately rather than paying a measurement
cadence it cannot use.

## End-to-end release gate H

Only run this after gates A-G pass; otherwise an apparently fast run may simply
be using a different budget.

Use the five ready modalities from `plan-adaptive-scheduling.md`: bulk SE, bulk
PE, scRNA, Flex, and scATAC, at `-t {8,32,64}` plus a small `-t {1,2,4}`
correctness set. Include dense-member, marker-window, sequential/inelastic,
plain, mixed regular/FIFO, warm-cache, and cold/source-bound cases.

For every modality and budget:

1. Establish an aggregate pinned oracle sweep from the same binary. Randomize
   run order; use at least three repetitions and five when within 3% of a gate.
2. Verify controlled-slot and CPU budget telemetry before accepting a timing.
3. Require exact input/output record counts and content-equivalent RAD output;
   retain byte-identical `-t 1` baselines and adversarial churn coverage.
4. Require median `auto` mapping time <=110% of the pinned oracle and no worse
   than the flat one-third split. A 95% uncertainty interval that crosses either
   boundary requires more repetitions, not a pass.
5. Require convergence within five moves and five seconds, and also within 20%
   of mapping time for short runs. A stable workload must have no resurvey after
   convergence.
6. Require total CPU seconds and peak RSS no worse than the best fixed split by
   more than 3% and 2%, respectively, unless the wall-time improvement is
   explicitly judged worth a documented trade.

Any one workload more than 10% off a correctly budgeted oracle falsifies the
current release claim. Do not average failures away across modalities.

## Recommended implementation sequence

### Milestone 0 — Make the budget and measurement experiment reproducible

- Add the authoritative `ExecutionPlan` and fix serial/no-feature, cgroup,
  scATAC, environment, and `-t 1` behavior.
- Add an aggregate pinned control and raw JSON telemetry.
- Complete gate A before collecting more performance numbers.

### Milestone 1 — Make actuation safe and observable

- Change resize traits and broker lifecycle methods to return errors.
- Add dynamic consumer floor and actual active-slot observations.
- Implement two-phase shrink acknowledgement before growth.
- Add failure injection and real multi-share churn tests.
- Complete gates B and C.

### Milestone 2 — Repair measurement before tuning

- Add upstream cumulative decode busy-time/component counters.
- Measure excluded consumer components and buffer occupancy.
- Replace fixed blackout completion with acknowledgement plus clean evidence.
- Complete gates E and G.

### Milestone 3 — Make the controller robust to wrong/noisy/nonlinear models

- Measure allocation-dependent response curves.
- Add epoch-scoped tabu with a nearest-safe fallback.
- Add robust ratification and uncertainty-derived hysteresis.
- Extend the discrete-event fake.
- Complete gates D and F plus the TB-08/TB-09 statistical gates.

### Milestone 4 — Validate the product, not just the crate

- Run correctness, churn, budget, oracle, convergence, cold-cache, FIFO, and
  overhead matrices.
- Complete gate H and the short-run/overhead gate.
- Treat the existing Flex `-t 32/64` measurements as two cells in this matrix,
  not as general validation.

### Milestone 5 — Documentation and release readiness

- Rewrite stale probe/supervisor docs and all override semantics.
- Make rustdoc pass with warnings denied in both feature configurations.
- Land or replace the two branch dependencies.
- Publish a release note for the changed `parallel=N` budget behavior only after
  its final semantics are settled.

## Verification performed during this audit

Passed:

- `cargo test -p thread-broker` — 16 controller tests and one doctest.
- `cargo clippy -p thread-broker --all-targets -- -D warnings`.
- `cargo clippy --all-targets --features rapidgzip -- -D warnings`.
- `cargo test --features rapidgzip` — 214 library tests passed, one ignored;
  fixture-dependent integration tests remained ignored in this command.
- `cargo test --release --features rapidgzip --test broker_adapters -- --ignored
  --nocapture` — all five real-adapter tests passed against the local fixtures.
- `RUSTDOCFLAGS='-D warnings' cargo doc -p thread-broker --no-deps`.

Failed as expected from stale documentation:

- `RUSTDOCFLAGS='-D warnings' cargo doc -p piscem-rs --no-deps --features
  rapidgzip` — five broken/private intra-doc links.

These passing tests do not contradict the findings above: the missing cases are
precisely budget propagation, asynchronous/failed actuation, multi-reader share
floors, producer measurement accuracy, buffer transients, nonlinear scaling,
and end-to-end modality coverage.
