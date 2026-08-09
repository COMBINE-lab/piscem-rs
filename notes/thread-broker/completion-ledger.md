# Thread-broker improvement completion ledger

This is the live execution record for `audit.md`. It is kept on
disk and deliberately ignored by git. Update it when work is started, completed,
measured, blocked, or invalidated; do not use it as a substitute for tests or
checked-in user documentation.

Started: 2026-08-06
Last updated: 2026-08-08
Baseline: `ac451f7` on `feat/thread-broker`

## Status legend

- `[ ]` not started
- `[~]` in progress
- `[x]` complete and its stated gate passed
- `[!]` blocked or failed; explanation required in the evidence log

## Milestone 0 — Authoritative budget and reproducible measurement

- [x] Define one immutable execution plan used by all mapping CLIs.
- [x] Fix serial and builds without `rapidgzip` to allocate the full effective
  budget to mapping.
- [x] Propagate `available_parallelism` to pool maxima and broker budget.
- [x] Remove scATAC's accidental rebinding of the total budget.
- [x] Define non-panicking `-t 1` behavior.
- [x] Define and test pinned/environment override semantics.
- [x] Add an aggregate fixed-producer control suitable for oracle sweeps.
- [x] Emit the effective execution plan in machine-readable telemetry.
- [x] Pass audit gate A.

## Milestone 1 — Safe, observable actuation

- [x] Make consumer and producer resize operations fallible/acknowledged.
- [x] Surface broker thread-spawn and thread-panic errors.
- [x] Respect the consumer floor in the initial split.
- [x] Model the real dynamic consumer floor from active reader groups.
- [x] Implement shrink acknowledgement before growing the other stage.
- [x] Add failure-injection and real multi-share churn tests.
- [x] Pass audit gates B and C.

## Milestone 2 — Reliable measurements

- [x] Make producer busy-time measurement accurate enough for control. The
  production adapter reads rapidgzip's feature-gated exact cumulative
  executing-region counter, so piscem starts no recurring sampler. The final
  monotonic snapshot protocol passed its direct 30-pair <=1% hot-path gate;
  feature-off code contains no clock read, counter update, or branch. Optional
  upstream thread-lifetime CPU accounting independently closes component
  reconciliation.
- [x] Account for or explicitly exclude auxiliary decode work.
- [x] Measure consumer state-construction and output-flush components.
- [x] Add clean-window/flow-transient evidence.
- [x] End blackout using acknowledgements plus clean evidence.
- [x] Pass audit gates E and G. Gate E passes with rapidgzip's exact
  cumulative executing-region signal across the full real input-archetype
  matrix. Both accounting-overhead gates and duration-based cap sub-gates pass.
  The formal persistent-source and temporary-lull Gate G halves pass at 25 and
  100 ms in all six cells apiece with canonical output equality.

## Milestone 3 — Robust controller and model

- [x] Replace permanent absolute tabu with epoch-scoped evidence and fallback.
- [~] Make ratification robust and distinguish inconclusive evidence. Interval
  comparison and explicit inconclusive outcomes pass the 1,000-seed synthetic
  acceptance thresholds. A pooled replay of the final bursty workloads shows
  that the proposed 95%-power real-noise gate is not statistically identifiable
  at the production horizon; see the final evidence entry.
- [x] Derive hysteresis from uncertainty across small and large budgets. All
  seven budgets pass the 1,000-seed stability and reopening gate.
- [~] Characterize and handle allocation-dependent stage scaling. The nonlinear
  simulation gate passes, and real short-run evidence added a bounded
  producer-pressure shrink veto plus a small-budget warm start. The full real
  five-split response curves remain.
- [~] Expand the fake to delayed/refused resize, queues, noise, lulls, and
  nonlinear rates without unsafe mutation. Failure injection, noise, lulls,
  nonlinear rates, and counter resets are covered.
- [x] Improve bounded, machine-readable reporting and terminal-state semantics.
- [~] Pass audit gates D and F and the TB-08/TB-09 statistical gates. Gate D,
  TB-09, the simulated portion of Gate F, and selected real normal-tier oracle
  cells pass. TB-08's synthetic gate passes, but its original pooled-real 95%
  power criterion is an information limit rather than a controller failure.

## Milestone 4 — Product-level validation

- [~] Correctness and adversarial churn across every modality. Generated real-
  pool churn passes; full end-to-end modality content comparison remains.
- [~] Budget telemetry validation before accepting any timing. The oracle runner
  rejects adaptive cells if any requested or observed trajectory, final split,
  terminal error, or peak occupancy violates the effective budget. Fixed cells
  validate their immutable pool configuration; independent runtime peak
  telemetry for fixed controls remains part of the overhead harness.
- [~] Oracle sweeps for bulk SE, bulk PE, scRNA, Flex, and scATAC at the planned
  thread counts.
- [~] Cold-cache, FIFO, mixed-input, source-bound, and inelastic coverage.
  FIFO-only, mixed regular/FIFO, and split-group cells now pass; cold-cache and
  real source-bound/inelastic cells remain.
- [~] Convergence, CPU, RSS, overhead, and short-run gates. CPU/RSS are captured,
  two opposing ~0.25-0.5 s bulk cells pass the short-run bound, and terminal
  states are honest. A builder-selectable freeze policy and configurable
  responsive probe cadence are implemented, with a local-only whole-broker
  fixed/default-responsive/sparse-responsive/freeze runner. The formal duration/input-count policy matrix
  remains.
- [x] Pass the formal long scATAC 64-versus-256 completed-progress publication
  gate with canonical equality and <=1% one-sided upper wall/CPU bounds.
- [~] Pass audit gate H. The sampler no longer blocks it, but the full
  correctness/performance and real Gates E-G matrices must run rather than be
  inferred from the preliminary cells.

## Milestone 5 — Documentation and release readiness

- [x] Rewrite stale probe/supervisor documentation and CLI help.
- [x] Make rustdoc pass with warnings denied in both feature configurations.
- [x] Settle and document `parallel=N`, environment, and aggregate-pin semantics.
- [!] Resolve or explicitly block release on the two unreleased dependencies.
  Both are pinned to reviewed commits, but neither has an acceptable released
  replacement; release remains explicitly blocked.
- [x] Prepare the behavior-change release note.

## Finding disposition matrix

This table is the cross-check that every audit finding has an implementation
path and that incomplete evidence is not hidden inside a completed milestone.

| Finding | Disposition | Remaining acceptance work |
| --- | --- | --- |
| TB-01 budget authority | Implemented; Gate A passed | None |
| TB-02 actual budget enforcement | Implemented; Gate B passed | Full Gate H repeats the observation across modalities |
| TB-03 surfaced actuation failures | Implemented; Gate C passed | None |
| TB-04 permanent tabu | Implemented; Gate D passed | None |
| TB-05 producer measurement | Exact upstream event-time counter adopted; direct no-measurable-overhead gate and the full real Gate E matrix passed | None |
| TB-06 allocation-dependent scaling | Partially implemented | Run five or more pinned splits per modality/budget; if cost varies by >10%, fit/validate empirical curves or extend the bounded correction |
| TB-07 flow/cap transients | Implemented; real Gate G passed | None |
| TB-08 ratification confidence | Synthetic gate passed; pathological real replay retained as a failed applicability check | Calibrate/replay at least 1,000 traces from the stable, non-injected real Gate F/H workloads; do not use alternating/source-lull traces as a stationary ratification population |
| TB-09 scaling/report semantics | Implemented; statistical gate passed | None |
| TB-10 robustness tests | Counter-reset and FIFO coverage implemented | Carry the portable feature matrix into CI |
| TB-11 reproducibility/docs | Implemented | Land or release the two pinned upstream revisions; this remains a release blocker |
| TB-12 overhead/short jobs | Same-split sampler overhead, scATAC fine publication, and two short cells passed; whole-broker A/B policies and local runner implemented | Complete the `{0.25,0.5,1,2,5}` s short-run matrix and 30-block `{5 s,1 min,10 min}` fixed/responsive/freeze cost matrix |
| Gate H | Selected normal cells passed | Run the complete normal modality catalog before release; comprehensive remains opt-in by risk |

## Remaining execution order

1. **Gate E complete.** The sampled fallback failed sparse/bursty
   decision-window accuracy, so production now reads rapidgzip's feature-gated
   exact cumulative executing-region counter and creates no sampler thread.
   Feature-off rapidgzip still has no hot-path clock/counter/branch; feature-on
   cleared its separate <=1% equivalence gate. See the final Gate E evidence
   entry below.
2. **Close Gates G, F, and TB-08.** Run the alternating-compressibility/lull
   fixture, then collect at least five pinned response points per real cell.
   Replay the observed variance through the 1,000-seed ratification test. Use an
   allocation-aware curve only where held-out prediction or >10% cost variation
   proves the current bounded correction insufficient.
3. **Close correctness and short-run coverage.** Use the implemented canonical,
   order-independent RAD digest command for every manifest job; add FIFO,
   cold-cache, source-bound, mixed-input, and the remaining
   duration cells. Counts alone are not an acceptable correctness pass.
4. **Run Gate H.** Use the documented local oracle runner for all five modalities at
   `-t {8,32,64}` plus `-t {1,2,4}` correctness cells. Keep raw JSONL, use at
   least three repetitions (five near a boundary), and enforce every per-cell
   wall/CPU/RSS/convergence/content gate without averaging failures away.
5. **Close the whole-broker overhead gate.** On stable real workloads, run at
   least 30 counterbalanced blocks of oracle-pinned/no-broker, responsive, and
   sparse responsive, and freeze-after-convergence at approximately
   `{5 s,1 min,10 min}`. Require exact
   canonical output and final-split equality. Report paired absolute seconds and
   fractional mapping-wall/process-wall/process-CPU ratios; every one-sided 95%
   upper bound must be <=1% as a broad regression backstop. Do not interpret 1%
   of aggregate process CPU as cheap: at 64 saturated threads it permits 0.64
   core. Primary direct lifetime accounting must bound controller+sampler CPU to
   <=0.001 core for responsive modes and <=5 ms total for freeze. Freeze sample
   count must stay within 25% across durations and absolute overhead must remain
   fixed rather than scaling.
6. **Release dependencies.** Replace both git pins with released crates (or
   record an explicit release exception) before shipping.

## Evidence log

### 2026-08-06 — Baseline

- `cargo test -p thread-broker`: 16 controller tests and one doctest passed.
- `cargo clippy -p thread-broker --all-targets -- -D warnings`: passed.
- `cargo clippy --all-targets --features rapidgzip -- -D warnings`: passed.
- `cargo test --features rapidgzip`: 214 library tests passed; fixture tests
  remained ignored by the default command.
- Release-mode real adapter tests: all five passed with local fixtures.
- Strict thread-broker rustdoc passed.
- Strict piscem rustdoc with `rapidgzip` failed on five stale/private links.
- Audit verdict: retain the control-law direction, but do not release before
  gates A-C and do not accept performance results before reliable budget and
  measurement gates.

### 2026-08-06 — Milestone 0 / gate A

- Replaced the three mutable/ad-hoc CLI budgets with the shared serializable
  `ExecutionPlan`. The requested budget, cgroup/cpuset-capped effective budget,
  mapping target, aggregate decode target, allocation kind, pin source, and any
  pre-clamp fixed request remain distinguishable.
- Bulk, scRNA, and scATAC now derive decoder headroom, mapping pool maximum,
  broker budget, startup targets, log fields, and `map_info.json` telemetry from
  that final plan. Actual parallel decoder handles are reconciled before work
  begins, returning unused reservations to mapping.
- Defined fixed controls: `--decoder parallel=N` and legacy
  `PISCEM_RAPIDGZIP_THREADS=N` are per decoder-capable input;
  `PISCEM_DECODE_SLOTS=N` is aggregate. Fixed allocations disable adaptation,
  reject zero/invalid/conflicting environment values, and visibly clamp to
  preserve one mapping slot. A one-slot effective budget always runs serially.
- Table-driven pure-planner tests cover requested budgets `{1,2,8,32}`, usable
  budgets including `{1,2,6,64}`, input counts `{0,1,2,17}`, serial/adaptive and
  both fixed allocation forms. They prove the effective-budget equality and
  both floors. Dedicated tests cover cap propagation, handle reconciliation,
  environment validation, one-thread behavior, and stable JSON fields.
- Evidence: default and `rapidgzip` planner suites passed (10 tests in each
  build); map-info tests passed (2); `cargo check` passed; strict Clippy passed
  for all targets with and without `rapidgzip`; `git diff --check` passed.

### 2026-08-06 — Milestone 1 / gates B and C

- Resize traits now return `Result`; producer adapters expose occupied execution
  slots separately from the requested limit. `start` and `finish` return typed
  errors for refused actuation, drain timeout, sampler spawn failure, and sampler
  panic. Runtime errors carry the partial report rather than fabricating a
  default report.
- Moves are two phase. The broker requests only the shrinking side, samples real
  `total_live`/`busy_workers` until the target is acknowledged, and only then
  grows the other side. A two-second default timeout leaves the pools at a safe
  under-budget split and returns an error. Reports record actual occupancy and
  peak controlled slots.
- The initial producer ceiling now subtracts the configured consumer floor.
  Each CLI computes a conservative paraseq floor from the collection's actual
  reader-group batch geometry and passes it to the broker. This accounts for the
  one-worker floor of every concurrently active share.
- The controller suite now has 23 tests. Injected failures cover consumer
  shrink, producer growth, producer rollback shrink, consumer rollback growth,
  a shrink that never acknowledges, sampler panic, and an initial producer
  request of `budget - 1` with an eight-thread consumer floor. Every runtime
  failure is surfaced with its partial report and remains at or below budget.
- Removed the raw-pointer mutation from the regime-change test; mutable costs
  now live under the simulation mutex, eliminating its data race.
- Portable real-pool tests use generated input only. The paraseq test churns
  `{1,2,8,24}` logical input groups at 100 ms, continuously monitors occupancy,
  and matches fixed-run record counts and sequence totals. The release-mode
  `rapidgzip` test does the same with generated multi-member gzip files and both
  real pools, waits on real busy-worker acknowledgements in both directions,
  enforces the execution-slot budget, and verifies average process CPU remains
  within the budget plus the documented one-CPU auxiliary allowance.
- Evidence: `cargo test -p thread-broker` passed (23 tests plus doctest);
  `cargo test --test paraseq_resize_safety` passed; release
  `real_pool_resize_safety` with `rapidgzip` passed; the real adapter refusal
  round-trip passed; strict all-target Clippy passed with and without
  `rapidgzip`.

### 2026-08-06 — Milestones 2 and 3 in progress

- The first audited decode adapter owned a dedicated 1 ms `BusySampler`,
  independent of the 25-250 ms controller cadence, and reported its thread plus
  upstream coordinator/scanner threads as auxiliary occupancy. A deterministic
  burst trace met Gate E's <=3% whole-run and <=10% p95 three-window error
  bounds. This provisional cadence was superseded by the phase-aware design in
  the 2026-08-07 evidence below.
- Producer adapters expose aggregate buffered progress. Windows with material
  fill/drain drift are rejected, and neither blackout nor ratification advances
  until it sees usable, stable-flow evidence after resize acknowledgement.
  Producer caps now require corroborating wall-clock duration, expire by
  wall-clock duration, ignore a lone lull, and report whether slack, source
  pressure, or both established the ceiling.
- Consumer busy metering now starts before per-callback state setup. Setup and
  output/merge flushing are separately accumulated, logged, and serialized in
  `map_info.json`, so the <2% exclusion decision can be made from data rather
  than assumed.
- Rejections now retain epoch, target, rates, uncertainty, and model share.
  They suppress a failed target only inside the unchanged epoch; a statistically
  persistent model/throughput shift opens a new epoch and permits revalidation.
  A deterministic two-regime test rejects target 11 in epoch 0, later reaches it
  as the true optimum, and stays within the five-moves-per-epoch bound.
- Ratification uses interval-separated block estimates and records inconclusive
  comparisons without creating rejection evidence. Cost-share uncertainty
  scales the effective deadband and resurvey distance by budget.
- Reports now contain bounded change-only requested/actual trajectories,
  dropped-sample counts, final phase and drift, terminal error, actual and peak
  occupancy, auxiliary occupancy, rejection evidence, epoch, uncertainty-
  derived thresholds, and cap reason. `converged()` requires a clean steady
  terminal state. The execution plan, broker report, and consumer component
  measurements are machine-readable in `map_info.json`.
- Current evidence: `cargo test -p thread-broker` passed 27 tests plus the
  doctest; the three `map_info` serialization tests passed with `rapidgzip`;
  `cargo check` passed with and without `rapidgzip`; `git diff --check` passed.

### 2026-08-07 — Statistical gates, documentation, and first oracle cells

- Seeded ratification exposed an over-conservative interval comparison: adding
  two 95% widths detected only 837/1000 true 10% regressions. Replaced it with
  a pooled one-sided 95% comparison. The corrected test has <1% false confident
  rejection on equal-rate blocks and >=95% detection on 1,000 independent
  traces at 2% window noise.
- Uncertainty hysteresis now has one-thread small-budget floors plus persistent
  reopening. For each budget `{2,4,8,16,32,64,96}`, 1,000 seeded 30-second
  stable traces meet the >=99% no-false-resurvey gate; every feasible 10%
  allocation shift meets the >=99% reopening-within-two-seconds gate.
- The coupled pipeline fake now makes both producer and consumer throughput
  deliberately sublinear in allocation. Auto delivers at least 90% of the
  discrete pinned oracle and settles within five moves. This passes the
  deliberate-sublinear simulation clause, not the real pinned-response-curve
  portion of Gate F.
- Producer cap tests at 25 and 100 ms establish a source ceiling after 300 ms,
  record its reason, and expire a cap caused by a one-second lull after 800 ms
  of contrary evidence. The burst sampler now cycles all required controller
  cadences `{25,50,73,100,137,250}` ms and keeps solved share movement within
  two points while retaining the whole-run and p95 error gates.
- Reports have an explicit bounded-trace/terminal-phase test. Machine telemetry
  now includes numeric mapping-pipeline seconds, excluding index load and output
  backpatching, so the oracle does not parse human logs or accept the old index-
  load confound.
- Added `scripts/thread_broker_oracle.py`: randomized same-binary fixed/auto
  sweeps, at least three repetitions, five fixed splits where possible, exact
  binary/commit/dependency provenance, raw JSONL, process CPU/RSS, telemetry
  rejection, correctness checks, and 10%-oracle/flat-one-third summaries. Its
  dry-run matrix passed.
- First real release-mode short-run cells, both bulk SE at `-t 8` with five
  fixed splits and three randomized repetitions, had exact matching counts and
  clean budget telemetry. PBMC passed: auto/oracle = 1.026 and beat the flat
  split. The 0.28 s SRR21186103 subset failed: auto stayed in survey at its
  one-slot opening split and was 1.226x the two-slot oracle. This is a genuine
  TB-12 short-run failure, not marked as passed. A budget-specific 25 ms / 100 ms
  small-run cadence is being tested against both opposing fixtures.
- User docs now describe the live broker and exact fixed-control semantics; the
  retired probe/supervisor study is clearly collapsed as historical. All three
  mapping CLI help pages are contract-tested. Strict rustdoc passes for the
  broker and piscem with and without `rapidgzip`. `paraseq` and
  `rapidgzip-core` are pinned to `40ae74fe...` and `dee928f...`; their lack of
  released equivalents remains an explicit blocker. Rapidgzip's required
  changes are merged upstream; paraseq's remain on a feature branch. A checked-in behavior-
  change release note is ready.

### 2026-08-07 — Small-budget short-run correction

- A 25 ms / 100 ms cadence alone did not repair the failing SRR cell: it reached
  steady state at one decode slot with a measured producer share near 5%, while
  the pinned surface still showed a 22% throughput gain at two. This is direct
  evidence of allocation-dependent service cost, not a convergence failure.
- Effective budgets up to eight now open at two decode slots where possible.
  A `Starved` producer may veto a model-requested *shrink* because runnable work
  is demonstrably queued behind the current limit. Pressure never chooses a
  larger target, so this is a bounded correction rather than a return to the old
  starvation hill climber. The veto is counted in machine telemetry. A coupled
  regression test with 10x producer under-reporting proves a true two-slot
  optimum is retained.
- Rebuilt the release binary and reran both randomized 18-run matrices unchanged.
  The previously failing SRR subset now passes: median auto 0.248 s versus
  two-slot oracle 0.252 s (`auto/oracle = 0.984`), exact counts, clean budget
  telemetry, and faster than the flat split. Its 0.25 s run ends honestly in
  `survey`, not false convergence.
- The opposing mapping-heavy PBMC subset also passes: median auto 0.459 s versus
  one-slot oracle 0.443 s (`auto/oracle = 1.037`), exact counts, and faster than
  flat. It starts at two and safely shrinks once to one; the short job ends in
  blackout/ratification and is again not mislabeled converged.
- These two cells satisfy the <=5% short-run bound for the tested duration and
  budget. They do not complete TB-12's full `{0.25,0.5,1,2,5}` second matrix or
  Gate H's modality/thread-count matrix, which remain in progress.

### 2026-08-07 — Final in-repository verification and handoff boundary

- The complete broker suite passes: 6 crate unit tests, 27 controller tests,
  and the doctest. The default workspace suite passes 227 library tests (one
  ignored) plus all enabled integration/doc tests; the `rapidgzip` suite passes
  228 library tests (one ignored) plus all enabled integration/doc tests.
- Strict Clippy passes for every target with and without `rapidgzip`. A final
  run caught and repaired a `needless_range_loop` in the cadence-ground-truth
  test; that test was rerun and passed. `cargo fmt --all -- --check` and
  `git diff --check` pass.
- Strict rustdoc passes for `thread-broker` and for `piscem-rs` with and without
  `rapidgzip`. The dependency still emits Cargo's pre-existing future-
  incompatibility notice for `proc-macro-error2 v2.0.1`; it is not a rustdoc or
  Clippy failure.
- The release-mode generated real-pool churn test passes. Against the local
  real adapter fixtures, all four ignored rapidgzip behavior tests pass; the
  separately enabled limit round-trip also passes. Together these cover
  throttled/starved, satisfied/ahead, sequential/inelastic, bounded busy time,
  and accepted limit behavior.
- The oracle now validates every adaptive requested and actual trajectory plus
  final and peak occupancy against the effective budget. It refuses to grant a
  correctness pass without one canonical, order-independent digest for the
  cell; equal record counts alone are retained as useful preliminary evidence
  but no longer satisfy Gate H. Existing auto and pinned short-run artifacts
  pass the strengthened telemetry validator.
- `benchmark_results/`, the local oracle/overhead runners, generated-data stress
  tests, fixture utilities, design notes, and this completion ledger are ignored.
  The audit, implementation, existing tracked regression tests, user
  documentation, and release note remain visible for normal source review.
- No remaining item is silently unplanned. Subsequent analysis found that the
  existing per-handle executing-worker signal is the correct lock-free source;
  the phase-aware accuracy/overhead work and revised Gate E plan are recorded
  below.

### 2026-08-07 — Phase-aware measurement and no-hot-path-change decision

- Rejected pool-wide `DecoderPool::stats()` as the integrator source. Reading it
  locks the member registry, and pool permit occupancy can remain high across a
  nonblocking task handoff. `BusySampler` now sums each cloned
  `DecoderHandle::stats().busy_workers`, whose relaxed atomic snapshot follows
  rapidgzip's executing CPU-region boundary.
- Replaced fixed 1 ms polling with trapezoidal integration at a deterministic
  jittered cadence averaging 3 ms only during survey/ratification. Draining,
  blackout, and steady monitoring use 25 ms; resurvey returns to calibration.
  Telemetry records both sample counts, mode changes, final mode, observation
  nanoseconds, and configured intervals. A settled-mode unit test proves the
  poll-rate reduction.
- The adversarial exact 5 ms on/off trace passed ten consecutive six-second
  trials at nominal 3 ms for controller windows `{25,50,73,100,137,250}` ms:
  <=3% whole-run error, <=10% p95 three-window error, and <=2-point solved-share
  drift. Nominal 4 ms failed (16.6% p95), so 3 ms is a measured boundary rather
  than a convenience constant.
- Added `PISCEM_FIXED_DECODE_MEASUREMENT=off|monitoring|calibration`, an
  environment-only fixed-aggregate validation hook. It preserves the exact
  mapping/decode split and binary while changing only sampler activity. The
  local overhead runner counterbalances all six within-repetition mode
  orders, requires canonical RAD equality, and gates both mapping wall time and
  process CPU using a one-sided 95% bootstrap upper bound of 1%.
- Formal results passed with 30 repetitions for one input, 60 for two inputs,
  and 30 for a lengthened 24-input fixture. Calibration wall/CPU upper bounds
  were respectively `{0.165%,0.510%}`, `{-0.021%,0.130%}`, and
  `{0.763%,0.162%}` relative to off. Monitoring bounds were
  `{0.298%,0.513%}`, `{0.134%,0.129%}`, and `{0.233%,0.000%}`. Canonical output
  and counts matched in every run. Median observation work ranged from 17-71 us
  for calibration and 4-20 us for monitoring over 0.25-0.57 s mappings.
- The sampler itself required no rapidgzip change. This initial decision was
  later superseded only for whole-run component reconciliation, not for the
  controller signal; see the next evidence entry.

### 2026-08-07 — Optional rapidgzip lifetime CPU accounting

- Created `feat/thread-cpu-accounting` from the prior reviewed rapidgzip pin
  `e4c7451`, implemented the compile-time optional `cpu-accounting` feature, and
  pushed it upstream. It is now merged on upstream `main` at `dee928f` and
  piscem pins that exact commit.
- The feature reads the per-thread CPU clock once at worker/auxiliary thread
  registration and once at thread exit, then performs one relaxed cumulative
  update at exit. No timing call, atomic operation, or conditional was added to
  decode task begin/end or inflate loops. Feature-off builds compile all
  accounting out; their stats fields are `None`. Running threads are
  intentionally omitted until exit, so this is component accounting rather
  than an exact decision-window busy counter.
- The direct rapidgzip gate used 30 paired, counterbalanced runs of a ~2 GiB
  decoded input with eight workers. Feature-on/off median CPU ratio was
  `0.989778` with one-sided 95% upper bound `0.998356`; wall ratio was
  `0.962264`, upper bound `0.991525`. Both passed the <=1% no-measurable-
  overhead criterion.
- A second 30-pair piscem gate compared separately compiled feature-off/on
  binaries at an identical fixed split with canonical RAD equality. Wall ratio
  was `1.002089` (upper `1.003515`) and CPU ratio `1.005089` (upper
  `1.007666`). All outputs matched. Median consumer busy + flush + completed
  worker CPU + auxiliary CPU explained `99.17%` of total process CPU, with zero
  CPU-clock failures.
- The full rapidgzip feature suite, feature/default Clippy and formatting pass.
  Piscem `cargo check --features rapidgzip-cpu-accounting` passes from the exact
  GitHub pin. Raw direct results are retained in rapidgzip's ignored
  `benchmark-results/cpu-accounting`; application results are in
  `/tmp/thread-broker-cpu-feature-overhead`.

### 2026-08-07 — Paired short-run oracle and first real response curves

- Corrected the oracle budget validator: a shrink-before-grow move may
  intentionally leave slots idle while acknowledgement is pending. Requested
  and actual trajectories must be non-negative and at most the effective
  budget; only final configured targets must sum to it.
- The runner now executes auto plus every fixed split as a counterbalanced
  repetition block and bootstraps the median of paired time ratios. It can
  enable fixed-split component measurement and writes median producer worker,
  auxiliary, sampled-busy, consumer-busy, setup, and flush components for every
  point.
- Five paired repetitions at `-t 8`, with exact canonical output equality,
  passed both short cells. Bulk SE auto/oracle median ratio was `0.9979`, upper
  95% `1.0141`, and auto/flat upper `0.8525`. PBMC auto/oracle was `1.0333`,
  upper `1.0983`, and auto/flat upper `0.8836`. Both CPU and RSS gates passed;
  there were no command or telemetry failures.
- Five fixed points expose material allocation dependence. Median completed
  worker CPU spans 288-398 ms for bulk and 320-404 ms for PBMC while consumer
  busy time stays within roughly 2%. This is Gate F evidence that the bounded
  live-pressure correction must be validated across the remaining modalities;
  it is not evidence for pretending service cost is constant.

### 2026-08-07 — Steady-state cost policies and whole-broker gate

- Added `SteadyStatePolicy::{Responsive, FreezeAfterConvergence}` to the public
  builder. Both execute identical warm-up, calibration, safe resize, blackout,
  ratification, and first clean steady validation. Freeze then terminates the
  controller and releases the only owner of the producer sampler, making its
  recurring cost exactly zero while explicitly giving up later regime changes.
- Added independent `steady_probe_interval(Duration)` configuration for
  responsive mode. It takes effect only after convergence. The producer adapter
  scales monitoring toward four observations per probe, clamped to the existing
  25 ms measured floor; survey/ratification remain at the jittered 3 ms cadence.
  A condition-variable stop path keeps `finish()` prompt with a 30-second probe.
- Piscem same-binary hooks are
  `PISCEM_THREAD_BROKER_POLICY=responsive|freeze-after-convergence` and
  `PISCEM_THREAD_BROKER_PROBE_INTERVAL_MS`. Telemetry includes policy, effective
  probe interval, controller lifetime/sample count, and explicit convergence
  shutdown. Lifetime rapidgzip CPU counters are refreshed at pipeline shutdown
  even when the producer sampler was released at convergence.
- Added `scripts/thread_broker_policy_overhead.py` and an example manifest. It
  counterbalances pinned/no-broker, default responsive, sparse responsive, and freeze blocks; validates
  canonical output, counts, budgets, final split, policy, convergence, and
  component accounting; and reports paired absolute and fractional mapping
  wall/process wall/process CPU cost with one-sided 95% bootstrap bounds.
- Unit evidence: 32 controller tests plus doctest pass. Dedicated tests prove
  freeze releases its producer before application finish, responsive retains it,
  zero probe intervals are rejected, a 30-second probe is promptly interrupted,
  and the producer monitoring interval changes at runtime. The formal real-data
  duration gate is still in progress and no whole-broker overhead number is
  claimed yet.
- The first 30-block whole-broker cell used a mapping-heavy bulk run with a
  3.798 s pinned mapping median and a real two-to-one-slot convergence. All 120
  runs had identical canonical output and final allocation. Relative to the
  fixed/no-broker baseline, default responsive, 5 s sparse responsive, and
  freeze had paired median mapping deltas of `+46.7 ms` (`+1.204%`), `+22.6 ms`
  (`+0.593%`), and `+23.7 ms` (`+0.623%`). Their process-CPU median/upper-95
  ratios were respectively `1.00587/1.00865`, `1.00100/1.00598`, and
  `1.00100/1.00837`.
- The deliberately strict universal <=1% upper-bound gate failed: mapping-wall
  upper bounds were 1.319%, 1.239%, and 1.041%; sparse/freeze process-wall
  bounds also missed by 0.092/0.394 points. This failure is retained rather than
  rounded away. Telemetry isolates a fixed convergence term: freeze stopped
  after 0.576 s and 23 controller samples; sparse used the same 23 samples over
  the run, while default responsive used 152. The next duration cell tests the
  user's stated production criterion: fixed absolute cost should become a
  decreasing fraction on minute-scale work, while responsive recurring cost
  must pass the direct <=0.001-core gate as well as the <=1% process backstop.

### 2026-08-07 — Administrative CPU gate corrected

- Clarified that the existing 1% CPU comparison uses aggregate process CPU. On
  a fully utilized 64-thread, ten-minute job it would permit 384 CPU-seconds,
  or an average 0.64 core, so it is retained only as a broad regression
  backstop rather than evidence that broker administration is cheap.
- Added lifetime CPU accounting for the broker controller and decode sampler.
  Each auxiliary thread performs exactly two thread-CPU-clock reads, at entry
  and exit; there are no new periodic clock reads in either loop and no changes
  to the rapidgzip hot path.
- The producer adapter is explicitly stopped and joined before the final report,
  including at freeze convergence, so sampler CPU is complete rather than a
  partial live snapshot. Telemetry records both CPU values and accounting
  failures.
- Updated the policy-overhead runner to reject missing accounting, report
  administrative CPU as absolute nanoseconds, fraction of one core over mapping
  wall, and fraction of aggregate process CPU. New primary gates are a one-sided
  95% upper bound of <=0.001 core for responsive modes and <=5 ms total for
  freeze; the paired <=1% whole-process gates remain to catch indirect effects.
- `cargo test -p thread-broker` passes 8 unit tests, 32 integration tests, and
  the doctest after this instrumentation.
- The rebuilt release binary completed 30 counterbalanced blocks (120 canonical
  runs) on the 3.8 s bulk cell with no output, allocation, convergence, or CPU-
  accounting rejection. One-sided 95% upper administrative CPU was `0.988 ms`
  for default responsive, `0.361 ms` for 5 s sparse responsive, and `0.334 ms`
  total for freeze. Expressed as fractions of one core over mapping wall these
  were `0.000255`, `0.0000933`, and `0.0000859`; all pass their new direct gates
  by wide margins. Median administrative CPU was only `0.00308%`, `0.00114%`,
  and `0.000996%` of aggregate process CPU.
- The direct pass does not waive the indirect regression backstop. All adaptive
  policies retained a fixed `36--40 ms` median mapping penalty, consistent with
  opening-split convergence rather than the sub-millisecond controller work.
  On this very short cell the one-sided mapping-wall upper bounds were `1.663%`,
  `1.558%`, and `1.744%`; process CPU passed for responsive/sparse but freeze
  narrowly missed at `1.0365%`. The overall cell therefore remains failed until
  the duration series demonstrates dilution on production-scale runs (or the
  fixed convergence term is reduced).
- Raw evidence: `/tmp/thread-broker-policy-bulk-4.6s-formal-cpu-v2`.

### 2026-08-07 — scATAC progress resolution and bounded nonlinear fallback

- The five-modality `-t 8` oracle matrix exposed a real nonlinear failure:
  scATAC auto collapsed to one decode slot and mapped 200,000 records in about
  16 s, while the pinned six-slot oracle took about 3 s. Producer CPU was only
  about 0.1 s; the cause was negative mapping-thread scaling, which a one-point
  isolated service-cost ratio cannot identify.
- Responsive mode now offers an opt-in probe of geometrically smaller consumer
  allocations after the model reaches the producer floor; piscem enables it
  only for scATAC, so other model-floor jobs pay no exploratory cost. Each point must
  confidently improve throughput, retained points compare against the previous
  measured point, and the first non-improvement reverts. Freeze is deliberately
  model-only and skips this cost. The deterministic negative-scaling test
  reaches producer 7 in three probes while freeze performs zero probes.
- Lowered only the scATAC paraseq reader batch from 16,384 to 1,024 records.
  This made consumer-shrink acknowledgement complete before the existing 2 s
  timeout. A fixed six-slot run improved from the prior 3.116 s median to 2.851
  s, so the liveness fix introduced no observed regression.
- Corrected a rate-estimation bias: zero-progress windows remain invalid for
  per-item service cost but now stay in throughput histories and ratification.
  Dropping them had conditioned slow, bursty scATAC estimates on observing a
  batch completion. A dedicated unit test fixes the zero/positive/zero mean.
- Replaced the temporary one-second/12-window workaround with the right
  consumer-layer fix. scATAC publishes completed-record progress every 64
  records while a
  responsive decision remains open, loads the cadence once per paraseq callback,
  returns to 256 in steady state, and uses drop-only publication for new batches
  after freeze convergence. A same-binary environment hook permits fixed-split
  A/B validation without adding a per-record configuration load.
- Decoupled fine progress from busy-time measurement after an exact-code A/B
  showed that publishing the full `(items, busy time)` pair every 32 records did
  not clear the <=1% overhead gate. The item counter can now publish at 64 while
  `Instant::now()` and busy-time atomics remain at 256. Unit tests prove that
  fine progress does not force fine clock reads and that freeze's maximum
  cadences publish only on batch drop.
- A real resolution sweep selected 64 rather than an arbitrarily fine cadence.
  All five 64-record runs followed 4→6→7→6 in three moves with no resurvey.
  At 128 records, two of five runs reopened, took six or seven moves, and made
  four or five nonlinear probes. Raw evidence:
  `/tmp/tb-scatac-progress{64,128}-auto-{1..5}`.
- After the publication fix, genuine short-range scATAC record-cost variance
  still required a bounded 12 × 25 ms nonlinear probe at `-t <= 8`; ordinary
  model moves remain four windows and larger budgets retain their 100 ms
  defaults. This is startup-only calibration, not recurring steady sampling.
  scATAC now opens adaptive `-t 8` at the measured safe midpoint (producer 4),
  retains producer 6, and rejects 7. Five independent 200,000-record runs all
  followed 4→6→7→6 in three moves with one revert, no resurvey, and
  1.78--1.95 s convergence; mapping times were 3.4708--3.5951 s. This is a
  material improvement over the earlier 1→4→6 path, but the short fixture
  still misses TB-12's <=5% safe-opening bound and is not marked complete.
- A provisional five-block, counterbalanced 2,000,000-record comparison against
  pinned producer 6 put adaptive/pinned paired median ratios at `1.050135`
  mapping wall, `1.049628` process wall, `1.049111` aggregate CPU, and
  `0.993242` peak RSS. The corresponding five-block bootstrap upper bounds were
  `1.100613`, `1.098583`, `1.102341`, and `0.995085`. Every adaptive run ended
  at producer 6 in three moves with no resurvey, but this approximately
  17-second cell is still too short to pass the 3% CPU gate. It is evidence for
  the remaining duration/warm-start work, not a Gate H pass. Raw evidence:
  `/tmp/tb-scatac-10x-paired-{auto,pin6}-{1..5}`.
- Added `scripts/scatac_progress_overhead.py`. Its initial randomize/reverse
  ordering was not actually counterbalanced: one 30-pair run placed fine first
  10 times and second 20 times. The runner now constructs exact 15/15 first
  positions, reports the counts, and rejects imbalanced evidence. Consequently,
  the earlier apparent pass is invalid and has been removed from the design and
  user documentation.
- The balanced 30-pair split-counter 32-versus-256 run retained canonical
  equality and produced median fine/generic ratios of `1.005160` mapping wall,
  `1.001493` process wall, and `1.001531` CPU. Upper-95 ratios were `1.020202`,
  `1.018018`, and `1.021696`; the three-second fixture is too noisy to satisfy
  the <=1% equivalence gate. The selected 64-versus-256 gate therefore remains
  pending a longer cell or enough additional balanced repetitions to resolve a
  sub-percent effect. Raw balanced 32-record evidence:
  `/tmp/scatac-progress-overhead-balanced-split-counters-32-vs-256`.
- Post-change verification passes: `cargo test -p thread-broker --all-targets`
  (10 unit and 36 integration tests), `cargo check -p piscem-rs` with and
  without `rapidgzip`, `cargo fmt --all -- --check`, and `git diff --check`.

### 2026-08-07 — 2-million-record scATAC split and policy assessment

- Built a 2,000,000-record, three-input concatenated-gzip fixture from the real
  scATAC input. Pinned and adaptive cells take roughly 13--32 seconds, so this
  measures behavior beyond the startup-dominated three-second fixture. Every
  accepted cell produced 2,000,000 reads, 528,410 mapped reads, and canonical
  RAD digest `293401b75cfd9a9003e3c83b0fe69070b68246f7c3b5a3ae69296674f36bb909`.
- A randomized three-repetition sweep covered every producer split 1--7 plus
  adaptive at total budget 8. Median mapping/aggregate-CPU seconds were:

  | Producer / mapper | Mapping s | CPU s |
  |---:|---:|---:|
  | 1 / 7 | 21.224 | 149.33 |
  | 2 / 6 | 21.129 | 127.88 |
  | 3 / 5 | 20.853 | 105.42 |
  | 4 / 4 | 17.396 | 70.79 |
  | **5 / 3** | **13.662** | 42.24 |
  | 6 / 2 | 17.865 | **37.03** |
  | 7 / 1 | 31.571 | 32.94 |
  | adaptive | 16.873 | 35.12 |

  Producer 5 is the wall-time oracle and was stable across 13.181--14.327 s.
  CPU and wall objectives differ: giving mapping fewer workers lowers aggregate
  CPU while eventually starving wall throughput. Raw sweep:
  `/tmp/thread-broker-scatac-long-oracle-64`.
- Adaptive was deterministic but wrong: all three sweep runs used
  4→6→7→6, converged in 1.779--1.804 s and three moves with no resurvey,
  and ended 23.5% slower than the producer-5 oracle. The geometric sequence
  never measures producer 5, so this is a search-resolution failure rather than
  noisy convergence.
- A separate four-block crossover-balanced policy comparison used pinned
  producer 5 as baseline. Every policy occupied every block position once.
  Medians and paired median ratios versus the pin were:

  | Strategy | Final producer / mapper | Converge | Mapping s | Paired mapping ratio | CPU s | Paired CPU ratio | RSS MiB |
  |---|---:|---:|---:|---:|---:|---:|---:|
  | pinned oracle | 5 / 3 | n/a | 13.981 | 1.000 | 43.18 | 1.000 | 1426.6 |
  | responsive, 25 ms | 6 / 2 | 1.892 s | 16.899 | 1.217 | 35.15 | 0.820 | 1413.4 |
  | responsive, 5 s | 6 / 2 | 1.816 s | 16.871 | 1.208 | 35.14 | 0.814 | 1418.5 |
  | freeze/model-only | 1 / 7 | 0.464 s | 19.731 | 1.422 | 137.67 | 3.209 | 1469.1 |

  Mapping ranges were 13.631--14.533 s pinned, 16.528--17.377 s responsive,
  16.304--17.767 s sparse responsive, and 18.934--20.965 s freeze.
- Sparse and default responsive made the same three moves, two nonlinear probes,
  and zero resurveys. Sparse/default paired medians were `0.99827` mapping,
  `0.99851` process wall, `0.99974` process CPU, and `1.00306` RSS: no observed
  throughput cost from the 5 s steady probe interval in this stable workload.
- Direct controller+sampler CPU makes the administrative cost explicit:

  | Strategy | Controller samples | Monitoring samples | Admin CPU | Fraction of one core | Fraction of process CPU |
  |---|---:|---:|---:|---:|---:|
  | responsive, 25 ms | 674 | 602.5 | 4.965 ms | 0.02936% | 0.01410% |
  | responsive, 5 s | 79 | 13.5 | 1.408 ms | 0.00835% | 0.00401% |
  | freeze/model-only | 22.5 | 8.5 | 0.401 ms | 0.00197% | 0.000282% |

  Administrative CPU passes the <=0.001-core gate for every policy. The large
  total wall/CPU differences are allocation effects, not broker-thread cost.
- The four-policy comparison holds the 64-record progress publisher present in
  every arm, including the pin, so it does not measure that mapping-worker cost.
  A separate four-pair, exactly 2/2 position-balanced fixed-producer-5 long A/B
  compared progress 64 against 256. Canonical output and counts agreed. Paired
  median fine/generic ratios were `1.009226` mapping wall, `1.006017` process
  wall, and `1.009200` aggregate CPU (median deltas 0.130 s, 0.085 s, and
  0.400 CPU-s). With only four pairs, upper-95 ratios were `1.052219`,
  `1.013928`, and `1.033157`; this is an estimate near the 1% boundary, not a
  formal pass. Raw evidence:
  `/tmp/scatac-progress-overhead-long-2m-64-vs-256`.
  This fixed-split experiment deliberately holds 64 for the entire run because
  no controller exists to switch phases; it is a worst-case continuous-publisher
  estimate. Adaptive uses 64 only through roughly 1.8 s of calibration and then
  returns to 256, so its publisher cost should dilute with run length. It also
  identifies a simple fixed-path improvement: absent the explicit validation
  environment override, pinned/no-broker scATAC should start directly at 256.
  Gate that change with identical canonical output and a non-regression relative
  to the current fixed default.
- Freeze is not an acceptable low-overhead scATAC policy in its current
  model-only form. All four measured runs moved once to producer 1, skipped
  nonlinear probing, and used 3.21x the oracle CPU. An additional retained
  freeze run confirmed the canonical digest; its large generated RAD/count
  payload was removed after hashing, while compact timing/map-info remains at
  `/tmp/tb-scatac-long-freeze-canonical-64`.
- The policy runner now uses four-period crossover rotation so each strategy
  occupies each block position exactly once per four repetitions; it reports and
  rejects positional imbalance. The interrupted preliminary producer-6 run is
  not used. Raw policy evidence:
  `/tmp/thread-broker-scatac-long-policies-pin5-64`.
- **Follow-up 1 complete: bounded interior refinement.** When geometric probing
  ends with a non-improving point and the last two proven points have an
  unmeasured discrete neighbor, the broker queues that neighbor exactly once.
  Monotone curves retain their prior probe count. A synthetic peak that exposes
  only endpoints 4, 6, and 7 now measures and selects producer 5.
- Four long, position-balanced auto/pin-5 pairs on the real fixture all produced
  the canonical digest and exact counts. Every adaptive run followed
  4→6→7→5, selected producer 5 in three moves/three probes with two retained
  improvements, converged in 2.305--2.455 s, and had no resurvey. Adaptive
  mapping was 14.185--14.752 s versus 13.959--14.565 s pinned; the paired median
  mapping ratio was `1.012329` and the four-pair upper was `1.050090`, passing
  the <=1.05 median and <=1.10 upper gates. Paired median CPU and RSS ratios
  were `0.987723` and `0.990516`. One process-wall sample had unrelated
  post-mapping delay, so mapping time remains the split-search gate. Raw compact
  evidence: `/tmp/tb-scatac-refinement-{auto,fixed}-{0..3}`; generated RAD/count
  payloads were removed after hashing.
- **Follow-up 2 in progress: freeze after full calibration.** Added the distinct
  `FreezeAfterFullCalibration` public policy and
  `PISCEM_THREAD_BROKER_POLICY=freeze-after-full-calibration`. The original
  freeze remains explicitly model-only. Deterministic tests prove that the new
  policy completes nonlinear/interior refinement, reaches producer 5, reports
  convergence, and drops its producer adapter; the suite passes 10 unit and 38
  integration tests. The policy runner now accepts a manifest-selected subset
  and includes the new mode.
- The first four-pair long gate deliberately failed: every full-calibration
  freeze run measured 4→6→7→5 but rejected 5 and restored 6, mapping in
  16.851--18.798 s. Shutdown itself was correct and total controller+sampler CPU
  was 1.664--2.228 ms, below the 5 ms gate. A current-binary responsive control
  reproduced the same rollback, ruling out freeze shutdown as the cause.
- A debug trace isolated the bias: the candidate-5 window immediately after
  producer 7 measured 51,942 records/s, but was compared with an earlier warm
  producer-6 baseline of 59,819 records/s. The physical transition grew mapping
  from one to three threads while the statistical comparison assumed an
  adjacent warm baseline. The controller now restores the retained best split,
  completes a nonlinear blackout, and only then probes the interior neighbor.
  Synthetic coverage asserts that physical restore. The rebuilt real gate is
  pending; failed raw evidence is retained at
  `/tmp/thread-broker-full-calibration-freeze` and
  `/tmp/thread-broker-current-responsive-diagnostic`.
- The retained-baseline correction alone was necessary but insufficient: the
  second four-pair gate kept producer 5 in only one of four runs. The three
  failures did restore 7→6 before 6→5 and stayed within five moves, but the
  ordinary 300 ms candidate blackout still admitted consumer cold-start. Admin
  CPU remained 2.031--2.413 ms. A corrected debug run measured the intended
  6→5 comparison as 81,315 versus 57,263 records/s and retained 5, confirming
  that the steady response is unambiguous once warm. Interior probes now receive
  a one-time doubled nonlinear blackout (600 ms for this scATAC configuration),
  without lengthening their measurement window or any steady-state cadence.
  The deterministic suite still passes; the release rebuild and third real gate
  are in progress. Failed second-gate evidence is retained at
  `/tmp/thread-broker-full-calibration-freeze-v2`.
- The third gate showed that a 600 ms blackout was also the wrong abstraction:
  it again retained 5 in only one of four runs. All failures respected the
  7→6→5 physical sequence and five-move bound; admin CPU was 2.110--2.771 ms.
  The successful run mapped in 12.561 s versus fixed runs at 13.370--14.435 s,
  so the settled candidate is decisively useful but its ramp length varies.
  Interior refinement now uses up to two *additional resolved 300 ms horizons*
  only when its first horizon fails to establish improvement. This bounded
  multi-horizon decision reuses 64-record progress, is reported as
  `nonlinear_probe_extensions`, and does not alter any steady-state polling.
  Failed third-gate evidence is retained at
  `/tmp/thread-broker-full-calibration-freeze-v3`; the fourth rebuilt real gate
  is pending.
- The fourth gate retained producer 5 in three of four runs, a material
  improvement, but still failed. The three accepted runs mapped in
  13.086--14.254 s versus fixed 13.274--14.350 s; the failed run returned to 6
  and mapped in 17.317 s. Two accepted runs required zero extensions and one
  required two; the rejected run exhausted two. The implementation error was
  that extensions accumulated the cold first horizon into every later estimate.
  Confirmation horizons are now independent, so a settled horizon is judged on
  its own signal, with at most four startup-only extensions. One otherwise valid
  run also recorded 8.351 ms administrative CPU (6.579 ms controller,
  1.772 ms sampler), above the 5 ms gate; this is retained as a gate failure and
  must not be averaged away. Raw evidence:
  `/tmp/thread-broker-full-calibration-freeze-v4`.
- **Follow-up 2 complete.** With independent confirmation horizons, all four
  balanced full-calibration-freeze runs followed 4→6→7→6→5, selected producer
  5 in four moves/three probes with no resurvey, and converged in
  3.132--4.125 s. Extensions were `0,1,1,3`, so no run reached the cap. Every
  run stopped monitoring, emitted exactly one final monitoring sampler
  observation, produced 2,000,000/528,410 read counts, and matched canonical
  digest `293401b75cfd9a9003e3c83b0fe69070b68246f7c3b5a3ae69296674f36bb909`.
  Against four position-balanced producer-5 pins, paired mapping median/upper
  ratios were `1.016213/1.036947`, passing the <=1.05/<=1.10 split gates;
  process wall was `1.016328/1.036036` and aggregate CPU
  `0.965353/0.990166`. Controller+sampler CPU was 2.250--2.582 ms with a
  2.582 ms upper bound, passing <=5 ms. Raw evidence:
  `/tmp/thread-broker-full-calibration-freeze-v5`.
- Next ordered task: default stable scATAC to sparse responsive monitoring while
  retaining `PISCEM_THREAD_BROKER_PROBE_INTERVAL_MS=25` as the normal responsive
  override; gate both modes on the current release workload and telemetry.
- **Follow-up 3 in progress.** scATAC now constructs its broker with a 5 s
  responsive steady interval before applying the existing environment override.
  Startup calibration, nonlinear probing, ratification, and resurvey cadence are
  unchanged. A feature-gated unit test asserts the 5 s default and the existing
  parser tests cover explicit positive intervals and reject zero. `cargo check`
  passes with and without rapidgzip; the release rebuild and real default/25 ms
  A/B telemetry gate are in progress.
- **Follow-up 3 complete.** An eight-pair, exactly 4/4 position-balanced direct
  comparison used default sparse responsive as treatment and explicit 25 ms
  responsive as baseline. All 16 runs selected producer 5 and matched canonical
  output. Sparse/default paired median/upper ratios were
  `0.992141/1.021036` mapping, `0.992385/1.020841` process wall, and
  `0.994118/1.019607` aggregate CPU, passing the <=1.05 upper gate. Median
  controller samples fell 566.5→144.5, producer monitoring observations
  431→10, and administrative CPU 4.902→2.623 ms; both policies remained below
  0.001 core. Sparse reports serialized 5 s, normal reports 25 ms, and neither
  stopped monitoring before application finish. Raw evidence:
  `/tmp/thread-broker-scatac-sparse-default-vs-normal-direct`.
- The earlier three-arm four-pair diagnostic is retained at
  `/tmp/thread-broker-scatac-sparse-default-vs-normal`. It had near-unity direct
  sparse/normal medians but a wide interval from one 15.9 s sample; that is why
  the balanced direct gate was run rather than declaring success from medians.
- Next ordered task: run the formal 30-pair long fixed-producer-5 scATAC
  64-versus-256 progress-publication gate with exact content/count equality and
  <=1% upper bounds for mapping wall, process wall, and aggregate CPU.
- **Follow-up 4 first formal gate failed.** All 60 runs completed with exact
  15/15 first-position balance, canonical equality, and fixed producer 5.
  Fine-64 paired medians were `+0.8797%` mapping, `+0.9203%` process wall, and
  `+0.8321%` aggregate CPU, but one-sided upper bounds were respectively
  `2.8821%`, `2.7690%`, and `2.6837%`, failing every <=1% gate. Raw evidence:
  `/tmp/scatac-progress-overhead-formal-long-2m-64-vs-256-final`.
- The fine path used one shared relaxed item atomic for every mapping worker, so
  64-record publication multiplied cache-line contention. scATAC progress is
  now written to cache-padded per-processor counters and summed by the single
  controller reader at its sampling cadence. Busy time remains on its separate
  256-record clock/atomic path, no per-record atomic load or clock read was
  added, and non-scATAC consumers retain the original meter. Dedicated tests
  prove external timer counters and sharded totals; the suite passes 10 unit and
  39 integration tests plus the piscem shard test. Release pilot/formal reruns
  are pending.
- A four-pair, exactly 2/2 position-balanced sharded pilot preserved canonical
  output and shifted fine/generic medians to `1.000559` mapping, `0.995048`
  process wall, and `0.998695` aggregate CPU, versus roughly `1.009` before
  sharding. As designed, four pairs cannot satisfy the formal interval gate;
  raw pilot evidence is `/tmp/scatac-progress-overhead-sharded-pilot-4`. This is
  sufficient to justify the required 30-pair rerun, not to declare a pass.
- The 30-pair sharded 2-million-record rerun was exactly 15/15 balanced and
  content-correct. Fine/generic medians improved to `1.003737` mapping,
  `1.003152` process wall, and `1.003721` aggregate CPU. Upper bounds were
  `1.012289`, `1.012210`, and `1.011713`: a large improvement over the shared
  counter, but still a formal failure by 0.17--0.23 percentage points. Raw
  evidence: `/tmp/scatac-progress-overhead-formal-sharded-2m-64-vs-256`.
- Because the remaining miss is interval precision around a 0.3--0.4% median,
  not a large point estimate, the next gate uses each 2-million-record input
  twice in one invocation (4 million records, roughly 28 s). It retains 30
  pairs, the same binary/counterbalance/canonical checks, and the unchanged <=1%
  upper bounds.
- The 4-million-record gate was also exactly balanced/content-correct but
  failed: paired medians were `1.005032` mapping, `1.005174` process wall, and
  `1.004578` aggregate CPU, with naive upper bounds near `1.02`. Longer runs did
  not reduce host/order variability. Raw evidence:
  `/tmp/scatac-progress-overhead-formal-sharded-4m-64-vs-256`.
- The runner's original bootstrap also discarded its crossover design by
  resampling all 30 ratios together, allowing imbalanced position mixes despite
  the real 15/15 schedule. It now resamples within each position stratum and
  combines their log-median ratios geometrically. This is the primary reported
  upper bound; pointwise unstratified paired medians remain visible. Reanalysis
  still leaves the prior datasets failed, so this correction does not manufacture
  a pass.
- Each sharded writer still used `fetch_add`, though no other writer touches its
  cache-padded counter. The timer now has an explicit single-writer cumulative
  publisher: it loads once per callback and uses relaxed stores at the 64-record
  cadence. Busy time and generic consumers retain their shared additive meter.
  A test proves cumulative semantics across successive timers. The 10-unit/
  39-integration broker suite, shard unit test, and rapidgzip build check pass;
  rebuilt measurement remains pending.
- The cumulative-store four-pair pilot preserved canonical output and had
  fine/generic paired medians `1.002281` mapping, `1.002477` process wall, and
  `1.002684` aggregate CPU. Two observations per position stratum are too few
  for a stable design-aware interval, so this justifies rather than replaces the
  formal rerun. Raw pilot evidence:
  `/tmp/scatac-progress-overhead-store-pilot-4`.

### 2026-08-07 — Ordered follow-up 4 complete

- The formal cumulative-store run completed 30 pairs / 60 mappings on the
  2-million-record fixture, with exactly 15 fine-first and 15 generic-first
  observations. Every run produced 2,000,000 input records, 528,410 mapped
  records, and canonical digest
  `293401b75cfd9a9003e3c83b0fe69070b68246f7c3b5a3ae69296674f36bb909`.
- Fine-64 versus generic-256 paired medians were `0.999332` mapping wall,
  `0.998945` process wall, and `0.999654` aggregate CPU. The crossover-aware
  position-stratified point estimates were `0.999332`, `0.998945`, and
  `0.999525`; their one-sided 95% upper ratios were `1.004188`, `1.004337`, and
  `1.004406`. All three pass the formal <=1% gate.
- The decisive hot-path change is a cache-padded, single-writer cumulative
  progress shard per scATAC processor. Publication uses a relaxed store at 64
  records while a decision is open. The controller sums shards at its existing
  sampling cadence. Generic consumers retain the shared additive counter and
  busy-time measurement remains at 256 records, so no additional clock read was
  introduced.
- Raw evidence: `/tmp/scatac-progress-overhead-formal-store-2m-64-vs-256`.
  This completes all four ordered follow-ups: bounded interior refinement,
  freeze after full calibration, sparse responsive scATAC steady monitoring,
  and the formal fine-progress publication gate.

## Local-only validation assets

The following assets were used to establish the gates above but are deliberately
not tracked by git. Their paths are explicitly ignored so a source commit cannot
accidentally include machine-local harnesses or generated evidence.

| Local path | Purpose |
| --- | --- |
| `scripts/thread_broker_oracle.py` and `.example.json` | Counterbalanced adaptive/pinned split sweeps, telemetry validation, canonical-output checks, and paired summaries |
| `scripts/thread_broker_overhead.py` and `.example.json` | Same-binary fixed-split producer-sampler cadence overhead gates |
| `scripts/thread_broker_policy_overhead.py` and `.example.json` | Fixed/default-responsive/sparse/freeze policy comparisons and direct administrative-CPU gates |
| `scripts/thread_broker_fifo_matrix.py` | One-shot FIFO-only, mixed-group, and split-group end-to-end correctness/message/plan matrix |
| `scripts/rapidgzip_busy_time_overhead.py` | Paired feature-on/off hot-path equivalence gate for the upstream exact busy-time counter |
| `scripts/thread_broker_measurement_gate.py` and `.example.json` | Real-archetype Gate E trace runner and cadence/component analyzer |
| `scripts/thread_broker_gate_fixtures.py` | Deterministic dense-burst and stored/source-dominated FASTQ fixture generator for Gates E/G |
| `scripts/scatac_progress_overhead.py` | Exact-position-balanced 64-versus-256 scATAC progress publication gate with stratified bootstrap intervals |
| `scripts/scatac_policy_geometry.py` | Position-balanced serial/pinned/adaptive scATAC batch/progress wall, CPU, RSS, liveness, and canonical-content gates |
| `examples/canonical_rad.rs` | Canonical, order-independent RAD digest helper used by the local runners |
| `examples/fastq_prefix.rs` | Deterministic FASTQ-prefix fixture generator |
| `tests/cli_help.rs` | Local CLI contract validation in both feature configurations |
| `tests/paraseq_resize_safety.rs` | Generated-data real paraseq pool resize/churn stress test |
| `tests/real_pool_resize_safety.rs` | Generated multi-member gzip rapidgzip/paraseq occupancy, acknowledgement, and CPU-budget stress test |
| `benchmark_results/` | Compact raw JSONL oracle results and machine-local run directories |
| `/scratch3/rob/tmp/thread-broker-*` and related `/scratch3/rob/tmp` fixtures | Large generated inputs, RAD/count outputs, map-info telemetry, and new retained pass/failure artifacts; older entries below retain their historical `/tmp` paths |
| `design.md`, `plan-adaptive-scheduling.md`, and this file | Local design, historical planning, and completion records |

Tracked regression coverage remains in the existing crate/module test suites
and already-tracked integration tests. The local-only tests above are valuable
for release qualification but depend on generated fixtures, optional features,
or machine-level timing/occupancy assumptions and are not part of the portable
source commit.

### 2026-08-07 — Pre-commit verification

- `cargo fmt --all -- --check`: passed.
- `cargo test -p thread-broker --all-targets`: 10 unit and 39 integration tests
  passed.
- `cargo check -p piscem-rs` and `cargo check -p piscem-rs --features
  rapidgzip`: passed.
- `cargo test -p piscem-rs --lib`: 229 passed, one existing fixture-dependent
  index-build test ignored.
- `cargo test -p piscem-rs --features rapidgzip --lib`: 236 passed, one existing
  fixture-dependent index-build test ignored.
- Strict Clippy passed for every thread-broker target and for piscem library and
  binaries with rapidgzip.
- `cargo build --release --features rapidgzip`: passed and refreshed
  `target/release/piscem-rs`.
- Cargo's pre-existing future-incompatibility warning for
  `proc-macro-error2 v2.0.1` remains informational; no command emitted a project
  warning or failure.
- The relevant implementation, existing tracked regression tests, audit,
  release note, and user documentation were committed on `feat/thread-broker`
  as `53388a6` (`feat(io): complete gated thread-broker implementation`). The
  local-only asset inventory above remained ignored and unstaged. This commit
  has not yet been pushed.

## Review-2 remediation execution

Started: 2026-08-07
Review: `review-2.md` against `53388a6`

- [x] Make every production broker lifecycle failure advisory without losing
  failure telemetry or output finalization; add injected regression coverage.
- [x] Run FIFO-only, mixed regular/FIFO, and split-group end-to-end correctness
  cells.
- [x] Document the minimum adopter contract and `Starved` shrink-veto exception;
  make the local validation tree pass standard all-target linting.
- [x] Scope scATAC reader batches and progress publication to the policy that
  needs them, with serial/fixed/adaptive A/B gates.
- [ ] Close real Gates E, G, F, and TB-08.
- [ ] Adopt and run the dual absolute/fractional whole-broker overhead gate.
- [ ] Run final Gate H.

### 2026-08-07 — Review-2 start

- Preserving the reviewer's uncommitted four-source-file fix and its local
  `.gitignore` entry. The fix correctly prevents `RunningBroker::finish()` from
  unwinding past RAD backpatching, but the first rescan found two remaining
  lifecycle gaps: `DecodeProducer::new(...)?` can still abort production mapping
  if the advisory busy sampler cannot spawn, and a start/panic failure with no
  partial `BrokerReport` is visible only in logs rather than `map_info.json`.
- Phase 0 will replace the lossy `Option<BrokerReport>` boundary with explicit
  advisory lifecycle diagnostics. Configuration and explicitly requested
  validation instrumentation remain fatal; production tuning failures fall back
  safely and remain machine-readable.

### 2026-08-07 — Review-2 phases 0 and 1 complete

- Production broker startup and completion are now represented by an
  `AdvisoryBroker` boundary whose `finish()` cannot propagate an error. Busy
  sampler spawn, controller spawn/initial resize, runtime resize timeout or
  refusal, and controller panic retain a structured stage/message in
  `thread_broker_failure` while mapping and output finalization continue at the
  last valid under-budget split. Explicit configuration and validation-only
  instrumentation errors remain fatal before the mapped result is accepted.
- Three injected rapidgzip tests cover measurement startup failure, initial
  resize refusal, and controller panic. The targeted tests, map-info diagnostic
  serialization test, and rapidgzip build check pass.
- The local one-shot FIFO harness ran all six bulk cells: FIFO-only, mixed
  regular/FIFO groups, and a regular-R1/FIFO-R2 split group, each under `auto`
  and explicit `parallel`. Every treatment completed without a lingering or
  failed writer, named every stream, used INFO for `auto` and WARN for the
  overridden explicit request, matched regular-file counts, and produced
  byte-identical RAD output. FIFO-only reconciled to serial `8/0`; mixed and
  split retained adaptive `6/2`, proving regular gzip inputs still opened the
  parallel path. Evidence: `/tmp/thread-broker-fifo-matrix/results.json`.
- Added local-only `scripts/thread_broker_fifo_matrix.py` to the validation
  asset inventory; it remains explicitly ignored and untracked.

### 2026-08-07 — Review-2 phase 2 complete

- The crate-level adoption guide now separates the minimum contract from the
  36-field policy surface: a truthful shared budget and floors, monotonic
  non-blocking busy counters, fresh-enough progress, truthful shrink
  acknowledgement, directional pressure, and an application-level advisory
  failure policy are required. Sampling, smoothing, blackout, ratification,
  caps, deadband, resurvey, regression tolerance, and nonlinear probes stay at
  defaults until named measurement failures justify changes. It also defines a
  concrete qualification set for a second adopter.
- `Producer::pressure`, `ProducerPressure`, `Starved`, and cap-history comments
  now state the bounded exception precisely: `Starved` may retain an occupied
  current allocation when the model asks to shrink, but can neither select nor
  increase a target. Source-bound and inelastic evidence remain growth caps;
  the cost-share model remains the only sizing mechanism.
- The ignored canonical RAD helper now supplies a default-feature stub and
  gates parity-only imports/functions. During validation, default all-target
  Clippy exposed and caused correction of missing `rapidgzip` guards in the
  advisory report-refresh blocks and constructor. `cargo clippy --all-targets
  -- -D warnings` and the same command with `--features rapidgzip` now pass in
  the full local validation tree.

### 2026-08-08 — Review-2 phase 3 complete

- scATAC now chooses pipeline geometry only after decoder-handle reconciliation,
  so a regular but non-gzip input that becomes serial does not inherit an
  adaptive-only cadence. The final reader batch and progress cadence are
  serialized as `pipeline_tuning` in `map_info.json`; same-binary positive-
  integer environment overrides remain available for validation.
- Fine 64-record progress is now adaptive-only. Serial and aggregate/per-file
  pinned runs publish at the generic 256-record cadence. Unit tests cover all
  allocation variants and explicit overrides. A missing broker cannot strand a
  fixed run at fine cadence.
- Measurement did not support scoping the 1,024-record reader batch away from
  fixed runs. On the exact two-million-record gate, 1K versus 16K had scoped/
  coarse medians of `0.9839` mapping wall, `0.9954` process CPU, and `0.9389`
  peak RSS for serial, and `0.9586`, `0.9690`, and `0.9874` for pinned. The
  smaller batch therefore has an independent long-run benefit and remains the
  scATAC default for every policy; the adaptive liveness reason is documented
  but is no longer its only justification.
- Adaptive 1K/64 versus forced 16K/256 had median ratios `0.6893` mapping wall,
  `0.8593` process CPU, and `0.9740` RSS. Every one of the coarse adaptive
  controls hit the same two-second consumer-shrink timeout; every scoped run
  completed the broker cleanly. The advisory failure path nevertheless
  finalized correct output in the negative control.
- The exact final long matrix used four position-balanced pairs per policy (24
  runs), 2,000,000 records/run, and one release binary. Every plan/tuning cell
  matched its request and every RAD/count result was canonical-equal. Serial's
  geometry-screening criteria are median <=1.05 plus a broad <=1.15 upper bound
  and <=10%/128 MiB RSS; pinned retains <=1.05; adaptive retains <=1.10 and
  converges within either 5 s absolute or 20% fractional wall. The later formal
  30-block whole-broker overhead gate remains the strict <=1% backstop. Evidence:
  `/tmp/thread-broker-scatac-policy-geometry-long-exact`.
- Short 200,000-record evidence, retained at
  `/tmp/thread-broker-scatac-policy-geometry-final`, showed the opposite serial
  startup trade: 16K was about 10% faster but used 6.5% more RSS. Because target
  runs are minutes to tens of minutes, the final policy follows the long-run
  result; the short-run screening envelope is explicitly <=20% median/<=25%
  upper rather than hidden. Intermediate 8K and 12K candidates cleared a 5%
  RSS screen but were performance-neutral with noisy six-pair intervals and
  were rejected; evidence is retained under the corresponding
  `/tmp/thread-broker-scatac-policy-geometry-{8192,12288}-pilot` paths.

### 2026-08-07 — Review-2 Gate E complete: native cumulative producer signal

- Added upstream rapidgzip feature `busy-time-accounting` on
  `feat/exact-busy-time-validation`. It integrates every existing decoder
  executing region, includes intervals still active at snapshot time, and is
  compiled completely out when disabled. A stable double snapshot affects only
  readers of telemetry; the feature-on hot path performs one monotonic clock
  read plus one relaxed wrapping-balance update at each existing begin/end
  boundary.
- The direct upstream gate used 30 randomized, counterbalanced pairs. Every run
  decoded the same `1,677,721,604` bytes/member result for eight repeated
  iterations with eight workers. Feature-on/off median CPU ratio was `0.995162`
  with one-sided 95% upper `0.998601`; wall ratio was `0.999995` with upper
  `1.009554`. Both satisfy the <=1% no-measurable-overhead gate. Raw evidence:
  `/tmp/rapidgzip-busy-time-overhead-exact`.
- The feature branch passed the full default and feature rapidgzip-core suites,
  focused active-interval coverage, strict feature Clippy, and formatting. It
  was pushed, merged with the concurrent Rust-version-only upstream change, and
  is now upstream `main` at `7b943ba3c8a669248ce70f3f958662859cd3e7c0`.
  Piscem pins that exact immutable revision; the local Cargo patch was removed.
- The first real Gate E matrix retained the 3 ms jittered sampler as the
  candidate signal and proved it is not release-quality. Whole-run error was
  <=2.5% on dense/bursty/mixed inputs, but p95 exact three-controller-window
  error was roughly 16--21% dense/mixed and 62--70% bursty; the short stored
  path was worse. Polling still faster would conflict with the low-overhead
  direction and would not make a sampled reconstruction exact.
- Production `DecodeProducer` therefore uses the native cumulative counter and
  creates no recurring busy sampler when the counter is present. The sampled
  implementation remains only as an older-upstream compatibility fallback and
  as an explicit fixed-split validation control. Native diagnostics have zero
  samples/observation cost and no sampler CPU; broker/controller cadence can no
  longer alias producer busy time.
- The final Gate E matrix ran three repetitions each of dense-member,
  marker-window, explicit sequential, bursty multi-member, mixed plain/gzip,
  and tenfold stored/source-dominated inputs. Observed final decoder paths were
  `DenseMembers`, `MarkerWindow`, `Sequential`, and `Stored`. At controller
  cadences `{25,50,73,100,137,250}` ms and five deterministic phase offsets,
  every nonserial cell had 0% whole-run and p95 three-window signal error and
  effectively zero solved-share drift. Sequential correctly emitted no
  producer signal. Event-time plus auxiliary components were at least 102% of
  measured producer CPU (the upper difference is non-CPU executing-region wall
  time), and consumer setup/flush was below 1.3% in every cell, so it remains an
  explicitly measured but excluded sub-2% component. Evidence:
  `/tmp/thread-broker-gate-e-real-native-v2`.

### 2026-08-08 — Review-2 Gate G/TB-08 real pilots and resulting corrections

- Added a bounded real-flow trace to `BrokerReport`: trajectory timestamps,
  cap/reason changes, and decision-open throughput windows. Steady-state
  throughput is not retained, so this adds no recurring telemetry work after
  convergence. `PISCEM_THREAD_BROKER_SAMPLE_INTERVAL_MS` supplies a positive
  same-binary active-cadence control for the 25/100 ms gate cells.
- Added a validation-only decoded-reader fixture hook and local flow runner.
  `PISCEM_VALIDATION_READ_LULL=MIB:MS[,MS...]` injects each finite lull once and
  leaves a clean recovery tail; an optional path match selects one stream in a
  paired run. The hook preserves rapidgzip's real decoder path and is inactive
  except for one environment lookup at input open. Unit tests cover strict
  parsing, finite duration exhaustion, and byte preservation. Local assets and
  generated alternating/Flex fixtures remain ignored.
- The first 72-million-read alternating-compressibility pilot ran for
  30.0--32.4 s/cell, kept canonical output equal, and held clean-window share
  bias to 1.9--2.3 points. It also established that this workload has a genuine
  persistent source ceiling near one slot; requiring its cap to clear after an
  injected lull was therefore an invalid combination of the persistent-source
  and transient-recovery sub-gates. Evidence:
  `/tmp/thread-broker-gate-g-pilot-v3`.
- A 100-million-pair Flex recovery pilot provided a nontrivial later demand of
  two to three decode slots and exposed a real cadence defect. The controller
  represented resurvey persistence as eight samples: at 25 ms it reopened 11
  times, while at 100 ms it needed roughly 1.9 s after the last lull to request
  the later-needed split. Clean share bias still remained below 2.1 points and
  canonical output was equal. This failed pilot is retained at
  `/tmp/thread-broker-gate-g-flex-pilot-v1`.
- Resurvey persistence is now 800 ms of wall time rather than a sample count.
  Persistent drift already contains a full clean cost window, so reopening
  reuses that evidence and asks for one fresh survey confirmation. This keeps
  the >=99% 30-second false-resurvey gate at every tested budget while leaving
  enough of Gate G's one-second recovery horizon for the confirming sample and
  actuation. The 11 unit plus 39 integration broker tests and strict Clippy pass.
- The first TB-08 replay incorrectly sampled independent residuals pooled
  across alternating regimes and failed with only 37.8% power. The corrected
  runner replays contiguous pre/post blocks and centers each block on the known
  null or 10% regression mean, preserving real within-block variance and
  correlation without relabeling an actual regime shift as noise. Four windows
  still had only 43.5% power. Ten windows passed 1,000 real-noise traces with
  0% false equal-rate rejection, 97.7% regression detection, and 0% isolated-
  zero rejection over 181 real windows. The production default is now ten
  windows (one second at 100 ms), and an isolated zero is explicitly
  inconclusive so it cannot create tabu evidence. Pilot evidence:
  `/tmp/thread-broker-tb08-real-pilot-v3.json`.
- Formal Gate G is split honestly: the alternating fixture validates clean
  bias and a persistent source ceiling, while the Flex fixture validates
  temporary-lull recovery at 25/100 ms. The three-repetition manifests are
  `scripts/thread_broker_flow_gate.formal.json` and
  `scripts/thread_broker_flow_gate.recovery.formal.json`; they remain local and
  untracked. Formal reruns follow the current release rebuild.

### 2026-08-08 — Final Gate F/H and whole-broker matrices prepared

- `scripts/thread_broker_oracle.formal.json` defines five modalities times
  budgets `{8,32,64}`, five pinned response points plus auto, and five rotated
  repetitions. Every fixed cell uses exact sampler-free producer accounting,
  and every admitted cell requires canonical RAD equality plus budget,
  occupancy, CPU, RSS, and convergence telemetry. The manifest and runner are
  local-only; outputs will be retained under
  `/tmp/thread-broker-gates-f-h-final`.
- Fixed-split measurement now accepts `PISCEM_FIXED_DECODE_MEASUREMENT=native`.
  It snapshots the same exact production counter and lifetime CPU components
  without starting the compatibility sampler; the oracle runner rejects a
  fixed cell if any sampler samples/CPU appear.
- `scripts/thread_broker_policy_overhead.formal.json` defines 30
  counterbalanced blocks for pinned, responsive, 5-second sparse-responsive,
  and freeze policies on approximately 5-second, 1-minute, and 10-minute
  versions of the same bulk workload. The runner now bootstraps both absolute
  deltas and fractional ratios. Acceptance is <=50 ms upper wall delta, <=1%
  upper mapping/process-wall and aggregate-CPU ratios, <=0.001 direct
  controller core for responsive modes, <=5 ms total freeze controller CPU,
  canonical equality, identical settled split, and fixed freeze work across
  duration. Generated 180-million and 1.8-billion-read inputs live on scratch3
  and all assets/evidence remain untracked.

### 2026-08-08 — Formal Gate G complete; TB-08 population boundary corrected

- The release-built formal alternating/source-bound matrix passed all six
  repetition/cadence cells under `/tmp/thread-broker-gate-g-real-v3`.
  Exact whole-run producer-share bias between paired control/lull executions
  was 0.08--0.76 percentage points, the inferred persistent ceiling was one
  slot, controlled occupancy never exceeded the eight-slot budget, and every
  canonical digest matched.
- The independent Flex recovery matrix passed all six cells under
  `/tmp/thread-broker-gate-g-flex-recovery-final`. Exact paired producer-share
  bias was 0.10--1.96 percentage points. Every lull that began after its paired
  clean control had established the later-required allocation recovered with
  zero observed delay, controlled occupancy stayed within budget, and all
  canonical digests matched. A lull before the clean control itself has learned
  that later allocation is deliberately not called a recovery failure.
- The final Gate G analyzer uses exact complete-work busy counters for its
  paired share check, the actually held later allocation for recovery, and only
  attributes caps whose onset is temporally associated with an injected lull.
  It does not require a sampled transient marker when an entire lull can fit
  between two 100 ms observations. These corrections define the measured
  quantities more precisely; none relaxes the five-point, one-second, budget,
  or canonical-content acceptance bounds.
- A 1,000-trace TB-08 replay that pooled both formal Gate G artifacts failed:
  false equal-rate rejection remained 0% and isolated-zero rejection remained
  0%, but ten-window detection of a true 10% regression was only 42.3% over 332
  windows (CV 19.3%). This is an expected information limit for the deliberately
  nonstationary alternating/source-lull fixture: no honest one-second 95%
  comparison can manufacture 95% power when adjacent blocks carry that much
  real regime variance. The earlier 97.7% pilot claim remains valid only for
  its lower-variance population and is not the final gate. Final TB-08
  calibration will use the stable, non-injected five-modality Gate F/H runs;
  the pathological failure is retained at
  `/tmp/thread-broker-tb08-real-final.json` as an applicability boundary.

### 2026-08-08 — Sparse resurvey audit and scaled Gate F/H fixtures

- Preflight of the first Gate F/H manifest rejected its nominal “formal” input
  sizes. At 32--64 threads, Flex finished in roughly one second and ended in
  blackout; bulk PE could finish before one post-warm-up requested/actual
  occupancy point existed. The interrupted 24-cell pilot is retained at
  `/tmp/thread-broker-gates-f-h-short-pilot-v1` and is not timing evidence.
- Repeated-gzip fixtures on `/scratch3/rob/tmp` now give the high-budget cells
  representative post-convergence work: bulk SE 180 million reads (8.58 s,
  converged at 2.10 s in its 64-thread smoke), bulk PE 20 copies (3.88 s,
  converged at 2.20 s), scRNA 100 million pairs (5.07 s, converged at 2.10 s),
  and Flex 100 million pairs (4.88 s, converged at 2.30 s). The existing scATAC
  fixture already runs substantially longer. Smoke RAD outputs were removed
  after telemetry inspection; fixtures and formal outputs remain local-only on
  scratch3.
- A code review before the scaled matrix found that the first sparse steady
  observation credited its entire five-second sleep toward the 800 ms resurvey
  persistence threshold. One noisy observation could therefore reopen a
  supposedly low-overhead broker. The first drift observation now earns only
  one active-cadence window; a second consecutive sparse observation spans real
  elapsed time and may satisfy the duration guard. A regression test pins that
  distinction. The 12 unit and 39 integration broker tests pass, and the release
  binary was rebuilt. The first three cells of the pre-fix scaled matrix are
  retained at
  `/scratch3/rob/tmp/thread-broker-gates-f-h-pre-sparse-fix-pilot-v2`.
- The oracle runner now has a fail-safe `--resume` mode. It retains only cells
  produced by the identical release-binary hash, rewrites prior failures out of
  the admitted matrix, and retries missing/failed cells. Without `--resume`, it
  refuses to truncate an existing `runs.jsonl`.
- Tightening the formal Gate G analyzer exposed that its persistent-source
  clause had checked only eventual cap presence, not the required <=1 s
  identification time. The first rebuilt matrix therefore fails despite its
  aggregate `gate_passed` field: direct cap telemetry appeared at 5--12 s and
  100 ms cells ended in ratification. It is retained at
  `/scratch3/rob/tmp/thread-broker-gate-g-real-v4-source-delay-failure`; a
  seven-run partial recovery rerun is retained at
  `/scratch3/rob/tmp/thread-broker-gate-g-flex-recovery-v2-pre-cap-fix-pilot`.
- Root cause: busy-derived slack and direct `Inelastic`/`SourceBound` pressure
  shared one `Caps::observe` call gated by stable downstream buffer flow, and
  cap telemetry was emitted only after a complete model snapshot. Pressure is
  now accumulated independently whenever the producer is not resizing; only
  busy-derived slack remains flow-sensitive. Cap changes are recorded directly
  from `Caps`, so the formal runner can enforce time-to-first `K+1` cap. A unit
  test proves three 100 ms direct-pressure observations establish the cap even
  without a clean cost window. The suite now passes 13 unit + 39 integration
  tests; a new release rebuild and both Gate G halves follow.

### 2026-08-08 — Counter-reset robustness

- Auditing the last TB-10 item found that `Work::delta` used saturating
  subtraction. If an adapter reset either cumulative counter, the broker would
  silently read zero work until the new counter caught the old value and could
  stall or misprice the controller. The controller now detects a decrease in
  either field immediately and returns structured
  `WorkCountersRegressed { side, previous, observed }` telemetry with its
  partial report. At the application boundary this remains advisory: mapping
  continues at the last safe split and output finalization is unaffected.
- Deterministic consumer- and producer-reset tests both pass. The complete
  release-mode broker suite is now 13 unit + 41 integration tests plus the
  doctest. Because this is a final-binary change, the formal real gates are
  rerun rather than inheriting evidence from the preceding build.

### 2026-08-08 — Monotonic exact rapidgzip accounting and final Gate E

- The first 300-million-record Flex rerun exercised the new piscem regression
  guard and caught a real producer decrease: cumulative busy time moved from
  `1,207,951,477` ns to `627,202,443` ns while completed items still advanced.
  The prior rapidgzip snapshot read two independently changing atomics twice;
  a cross-counter ABA could therefore admit a mismatched live-worker snapshot,
  transiently overestimate the total, and regress at the next observation.
- A completed-region-only alternative was implemented and measured rather than
  accepted on intuition. It was monotonic and passed the direct overhead gate,
  but the independent high-resolution sampler diagnostic showed 25--100 ms p95
  lag errors of 17--313% on real sparse/bursty/stored paths. It was rejected as
  the wrong control signal. The diagnostic is retained at
  `/scratch3/rob/tmp/thread-broker-gate-e-real-monotonic-v3`; its passing direct
  overhead result at `/scratch3/rob/tmp/rapidgzip-monotonic-gate/overhead-final`
  does not qualify the rejected design.
- The final upstream design retains exact live-region accounting and adds one
  feature-gated packed transition state. Its low 32 bits count active boundary
  transitions and its high 32 bits form a generation. A snapshot is accepted
  only when the transition word is unchanged with no active transition and the
  existing busy/balance values are stable and plausible. This detects both
  partial updates and complete ABA cycles. With the feature disabled, all
  clocks, counters, transition operations, and related branches compile out.
- The full upstream feature-on and feature-off suites passed, including a
  four-worker, 10,000-region-per-worker concurrent monotonicity stress test.
  Strict Clippy is blocked only by 14 pre-existing `collapsible_if` findings
  from the newer toolchain. The direct hot-path gate used 30 counterbalanced
  pairs with identical stdout and passed: CPU ratio median `0.9957667`,
  one-sided 95% upper `0.9985809`; wall ratio median `0.9982332`, upper
  `1.0000000`. Raw evidence is retained at
  `/scratch3/rob/tmp/rapidgzip-epoch-gate/overhead-final`.
- Upstream branch `fix/monotonic-busy-time-accounting` was pushed at `85d90c6`,
  merged to rapidgzip `main`, and pushed as merge commit
  `276a41f77fb927e24cb0898a638a08b21eb048c6`. Piscem is pinned to that merge,
  and its final release binary was rebuilt from the remote revision.
- A final remote-pinned Gate E run passed all 18 cells (six archetypes by three
  repetitions) at
  `/scratch3/rob/tmp/thread-broker-gate-e-real-epoch-v5-final`. Sequential input
  produced a truthful zero signal; DenseMembers, MarkerWindow, Sequential,
  and Stored paths were all observed; producer-component coverage ranged from
  1.02 to 2.00; and excluded consumer work remained below 1.5%. The accounted
  signal's algebraic self-comparison is not treated as independent accuracy
  evidence: acceptance also rests on the rejected-signal diagnostic, upstream
  concurrency stress, direct A/B, CPU-component reconciliation, and sustained
  real-workload monotonicity.

### 2026-08-08 — Final Gate G and recovery-fixture correction

- The final remote-pinned persistent-source matrix passed all six cells at
  `/scratch3/rob/tmp/thread-broker-gate-g-real-v7-final`. Persistent caps were
  identified in `0.376--0.400` s, exact paired producer-share bias was
  `0.0003--0.0029` percentage points, controlled occupancy never exceeded the
  eight-slot budget, and every canonical digest matched.
- Recovery pilots v7--v10 were retained rather than overwritten. They exposed
  two experimental-design errors: schedules were decoded MiB rather than
  compressed GiB, and a lull was called a recovery event even when the
  perturbed run did not hold the clean control's required allocation at event
  onset. Compound early lulls could therefore delay initial convergence and be
  misreported as 2--21 s recovery. The 300-million-read R1 stream was measured
  exactly at `112,042,871,190` decoded bytes (`104.35 GiB`).
- The final analyzer requires an event to begin after the clean control reaches
  its retained allocation and while the perturbed run itself holds that
  allocation. It still fails a cell with no such event and retains the one-
  second recovery bound. A single isolated one-second lull at `95 GiB` leaves
  more than `9 GiB` for post-event observation and avoids contaminating the
  precondition with earlier validation events.
- The final v11 matrix at
  `/scratch3/rob/tmp/thread-broker-gate-g-flex-recovery-v11-final` passed all
  six cadence/repetition cells. Every event was eligible, allocation recovery
  and harmful-cap clearance had zero observed delay, exact paired share bias
  was `0.00024--0.0144` percentage points, occupancy stayed within eight slots,
  and all canonical outputs matched. Unchanged v8 controls were reused only
  after exact binary-hash, command-prefix, cadence, repetition, successful
  telemetry, and canonical-digest validation.

### 2026-08-08 — Validation tiers and normal-tier switch

- Exhaustive qualification proved too slow for an interactive update. Local
  manifests now name three intentional tiers; comprehensive evidence remains
  available but the remainder of this run uses normal:

  | tier | oracle F/H | policy overhead | intended use |
  | --- | ---: | ---: | --- |
  | light | 32 runs: bulk SE + Flex, t32, auto + three pins, four position-balanced repetitions | 16 runs: four policies x four position-balanced 5 s blocks | routine implementation feedback; common opposing shapes only |
  | normal | 120 runs: all five modalities, four t8/t64 adaptive runs and four position-balanced auto + three-pin t32 blocks | 84 runs: four policies x ten 5 s blocks, ten 1 min blocks, and one descriptive 10 min block | default pre-merge qualification and the tier used below |
  | comprehensive | 540 runs: all five modalities/budgets, five pins + auto, six position-balanced repetitions | 360 runs: four policies x 30 blocks at 5 s/1 min/10 min | opt-in release/nightly or controller/measurement changes |

- Normal retains statistical upper-bound gates for 5 s and 1 min policy cells.
  Its single 10 min block is explicitly labelled descriptive; it measures the
  absolute/fractional trend and direct administrative CPU but is not presented
  as a 95%-confidence equivalence result. The full 30-block long-duration gate
  remains comprehensive.
- Oracle crossover blocks must now cover every mode position. The earlier
  three-repetition/four-mode long follow-up was invalid because adaptive
  scATAC always followed pin-4; the manifests use four modes/four repetitions
  in light and normal, and six modes/six repetitions in comprehensive. Normal
  retains a <=10% oracle upper bound and adds a <=5% median bound relative to
  the simple one-third baseline; requiring exactly 0% against a near-tied
  baseline was below the resolution of the measured workloads.
- Pathological recovery, full Gate E archetypes, FIFO/mixed-input, cold-cache,
  full short-duration grids, and direct rapidgzip A/B are opt-in suites selected
  by risk. They remain mandatory when code touches their associated input,
  accounting, lifecycle, or controller behavior, rather than on every ordinary
  update.

#### Per-policy overhead table

The completed 14--20 second scATAC crossover measurements quantify every
shipping method today. They do not replace the unrun normal duration series;
that series remains the gate for making a duration-scaling claim.

| policy | measured mapping wall | measured aggregate CPU | direct administrative CPU | interpretation |
| --- | ---: | ---: | ---: | --- |
| fixed/no broker | 13.981 s; 1.000 | 43.18 s; 1.000 | no broker thread | pin-5 oracle baseline |
| responsive, 25 ms | 16.899 s; 1.217 paired ratio | 35.15 s; 0.820 | 4.965 ms; 0.02936% core | older search selected 6/2; low direct cost did not compensate for split error |
| sparse-responsive, 5 s | 16.871 s; 1.208 paired ratio | 35.14 s; 0.814 | 1.408 ms; 0.00835% core | same split, 72% lower administration; closest to original sparse-sampling design |
| model-only freeze | 19.731 s; 1.422 paired ratio | 137.67 s; 3.209 | 0.401 ms; 0.00197% core | near-free but chose the wrong 1/7 split |
| full-calibration freeze | 1.016 median / 1.037 upper versus oracle | 0.965 median / 0.990 upper | <=2.582 ms | retained calibration quality, then stopped recurring work |

The later direct eight-pair sparse/25-ms comparison held the controller and
selected split constant: sparse/normal mapping was `0.9921/1.0210`
median/upper-95, aggregate CPU `0.9941/1.0196`, and administrative CPU fell
`4.902 -> 2.623 ms`. Remaining normal-duration fields are intentionally not
fabricated: 5-second, 1-minute, and descriptive 10-minute cells remain opt-in
at `/scratch3/rob/tmp/thread-broker-policy-overhead-normal-final`.

### 2026-08-08 — Controller landing checkpoint

- Final focused verification passed: `cargo fmt --all -- --check`,
  `git diff --check`, strict all-target Clippy both without rapidgzip and with
  `rapidgzip-cpu-accounting`, 15 controller unit tests, 41 control-law
  integration tests, and 243 piscem library tests (one fixture test ignored).
  The feature-off check caught and fixed a misplaced scATAC helper cfg. The new
  regression test proves that an adjacent grant cannot confirm the provisional
  cap it inherited without headroom.
- The final remote-pinned, CPU-accounting release binary passed the balanced
  scATAC t32 normal cell at
  `/scratch3/rob/tmp/thread-broker-gates-f-h-long-normal-v5-final`: canonical
  equality, convergence, median adaptive wall `13.594 s` versus pin-4
  `13.776 s`, adaptive/oracle upper-95 `1.0418`, CPU ratio `0.9868`, RSS ratio
  `0.9995`, and about `0.1 ms` controller CPU.
- The final cap-confirmation fix passed the balanced scRNA t32 normal cell at
  `/scratch3/rob/tmp/thread-broker-gates-f-h-long-normal-v6-final`: canonical
  equality, zero resurveys, 0--1 moves, convergence in `0.3--2.1 s`, adaptive
  wall `20.555 s` versus pin-11 `19.902 s`, oracle upper-95 `1.0760`, median
  one-third-baseline ratio `1.0328`, CPU ratio `0.9770`, and RSS ratio `1.0156`.
- A final selected persistent-source Gate G normal rerun passed at
  `/scratch3/rob/tmp/thread-broker-gate-g-real-normal-final`. The unchanged
  temporary-lull path retains the complete six-cell v11 pass at
  `/scratch3/rob/tmp/thread-broker-gate-g-flex-recovery-v11-final`; a redundant
  final-binary recovery rerun was stopped while canonicalizing its multi-GiB RAD
  after the control completed because the cap change only makes confirmation
  stricter and the already-passed recovery behavior was untouched.
- TB-08 remains an honest statistical applicability boundary. Replaying the
  pooled final bursty traces gave 77.4% detection with 10 windows; increasing
  to 20 windows reduced detection to 3.8% as the pooled CV rose to 38.9%.
  False equal-rate and isolated-zero rejection both remained 0%. This is not
  evidence for increasing production sampling: pooling heterogeneous startup,
  move, and modality regimes violates the stationary comparison population.
  The controller retains the synthetic 1,000-seed gate and ten-window default;
  any future real gate must define a stable per-modality epoch first.
- Landing deliberately does **not** wait for the 84-run normal policy-duration
  matrix. Remaining runnable catalog: normal 5 s/1 min/descriptive-10 min policy
  duration series; the unexecuted normal oracle cells needed for a complete
  five-modality release claim; comprehensive 540-run oracle and 360-run policy
  matrices; cold-page-cache harness; complete short-job grid; and risk-triggered
  Gate E/G, FIFO/mixed-input, scATAC geometry/publication, and rapidgzip direct
  A/B reruns. Light is the routine controller loop, normal is pre-merge/release
  qualification, and comprehensive is opt-in nightly or risk-triggered work.

### 2026-08-09 — Reviewer scATAC t32 failure and corrective gate

- Reproduced the reviewer's exact adaptive `map-scatac -t 32` failure with the
  merged release binary at
  `/scratch3/rob/tmp/tb-scatac-t32-repro-agent-a`: 10.955 s mapping, opening and
  final producer 4, zero moves, and zero nonlinear probes. The immediate cause
  was more specific than the review's guard hypothesis: piscem configured
  `min_producer_slots = 4` and disabled nonlinear probing above t8. The solver's
  internal target was therefore clamped to four even though the diagnostic
  unconstrained ideal was one. The ordinary deadband then correctly settled at
  the configured floor; no controller escape path was eligible.
- Expanded the real low-end fixed surface before selecting a correction. On the
  unchanged binary, t32 producer 1 mapped in 9.545 s and producer 2 in 8.348 s;
  the reviewer's balanced producer-2 median was 9.216 s. At t64, producer
  1/2/4 mapped in 3.943/4.111/4.286 s. Producer two is therefore the t32 oracle,
  is only 4.3% behind producer one at t64, and improves on the old producer-four
  floor at both large budgets.
- Intermediate real-data trials rejected broader changes. Moving from the old
  opening to producer one found the right model direction but retained a
  same-split startup/resize penalty at t32; lengthening baseline horizons and
  globally changing progress resolution did not remove it and carried wider
  modality risk. Starting at producer two removes that resize/cold-pool effect
  without additional sampling. The final large-budget policy therefore opens
  at 30/2 for t32 and 62/2 for t64, sets the measured floor to two, retains the
  generic 256-record progress cadence and 5 s sparse steady probe, and keeps the
  optional nonlinear search disabled. The ordinary model remains free to grow
  above two if the workload changes.
- Corrected the independent generic one-way-door weakness for opt-in nonlinear
  policies. Their floor-directed model move is now measured: a conclusive keep
  completes the response experiment; a regression or inconclusive floor
  discriminator restores the opening and starts the existing bounded
  upward/local search. Ordinary model moves retain the established
  confident-regression-only rollback rule. Unit and synthetic control-law tests
  cover both a valid mapping-heavy floor move and the t8 negative-consumer-
  scaling rollback/search path.
- Final same-binary real gates used release binary SHA-256
  `8fcb72290578a20f5798a0f940efa07ee47dc9952f775d5a5aa367e8542d2762`, all
  on the 2,000,000-read fixture:

  | budget | final mapping wall | selected path | convergence | controller CPU | comparison |
  | ---: | ---: | --- | ---: | ---: | --- |
  | 8 | 13.254 s | 4 -> 1 -> 4 -> 6 -> 7 -> 6 -> 5; 3 response probes | 4.911 s | 0.796 ms | faster than the reviewer's 13.76 s producer-5 oracle median |
  | 32 | 8.702 s | opens/holds producer 2; zero moves/probes | 0.300 s | 0.085 ms | 5.6% faster than the 9.216 s producer-2 oracle median and 24.4% faster than the 11.51 s failed adaptive median |
  | 64 | 3.943 s | opens/holds producer 2; zero moves/probes | 0.300 s | 0.071 ms | equal to the measured 3.943 s producer-1 best run; producer 2 remained inside the 8% split gate |

  Evidence is in `/scratch3/rob/tmp/tb-scatac-t{8,32,64}-fixed-agent-e`.
  Every run consumed 2,000,000 reads, mapped 528,410, reported no terminal
  error, and produced canonical digest
  `293401b75cfd9a9003e3c83b0fe69070b68246f7c3b5a3ae69296674f36bb909`.
- Final verification passed `cargo fmt --all -- --check`, `git diff --check`,
  `cargo check --no-default-features`, strict all-target Clippy with and without
  `rapidgzip-cpu-accounting`, all 16 thread-broker unit tests, all 42 synthetic
  control-law tests, and 243 piscem library tests (one fixture test ignored).
- The review's secondary no-op-resurvey `time_to_converge` reset and sensitivity
  of a one-thread resurvey distance remain separately cataloged improvements;
  neither caused this failure, and neither was mixed into this scoped fix.

### 2026-08-09 — Policy-knob review and self-bracketing opening

- Accepted the central finding in `policy-knobs.md`: a measured
  optimum is an opening hint, not a safety constraint. Work was done on
  `feat/thread-broker-opening-bracket`, branched from `dev` with the reviewed
  t32 correction cherry-picked as `0338bdc`. scATAC now leaves the shared-pool
  producer safety floor at one at every budget, opens with four producer slots
  at t8/t32/t64, and contains no `budget <= 8` controller-policy branch.
- Replaced four public nonlinear-search configuration knobs with the opt-in
  `OpeningPolicy::{Fixed, Bracket(OpeningBracketConfig)}`. The generic bracket
  defaults to at most three measured candidates, 200 ms of evidence per point,
  and a four-second total decision budget. It begins only after a stable model
  answer disagrees with the opening, measures the model answer first, and if
  that loses restores the opening and checks adjacent candidates first away
  from, then toward, the rejected model. An optional point must be both
  statistically separated and more than 5% faster before it is retained.
- The bracket is startup-only and is never rearmed by resurvey. Responsive and
  full-calibration freeze may run it; model-only freeze explicitly skips it.
  Model/opening agreement records `model_agreed` with zero points, samples and
  wall time. Public telemetry reports `points_measured`, `samples`,
  `wall_nanos`, and a bounded outcome. The ordinary resize blackout remains
  separate from the evidence horizon; reusing the horizon as a blackout had
  added roughly 400 ms to the measured t8 path without adding evidence.
- Candidate generation now respects both configured safety floors. This was
  caught during the final generic audit: hard-coding a candidate floor of one
  would have let a future adopter with a legitimate higher producer floor probe
  below its contract, even though scATAC itself could not expose that error.
- The paired synthetic sanity gate uses the same opening (producer four) and the
  same allocation-independent model answer (producer one) for opposite response
  shapes. The mapping-heavy t32-like shape retains one; the negative-scaling
  t8-like shape rejects one and proves producer five. Separate tests cover zero-
  cost model agreement, point/sample/wall bounds, deadline exhaustion, freeze
  semantics, no bracket rearm on resurvey, invalid configuration, and safety-
  floor-preserving candidate generation. Generic bulk/scRNA/Flex broker config
  is asserted to retain `OpeningPolicy::Fixed`, so those modalities cannot pay
  bracket work.
- The local, intentionally untracked oracle runner now refuses to call a fixed
  winner an oracle unless the swept pin set contains a measured point on both
  sides. Boundary winners are reported with `oracle_mode: null`, the candidate
  is retained for diagnosis, and the cell fails qualification as unbracketed.
  Three local unit cases cover an interior winner and both boundaries. The
  runner, tests, manifests, fixture data, and benchmark outputs remain local
  validation infrastructure and are documented here rather than tracked.
- Final real validation used release binary SHA-256
  `481b2be0163aa15c9feed1e1ac4f8c31d9f0c309447cf70809e3fa78f7380126`
  on the two-million-read / 528,410-mapped scATAC fixture:

  | budget | mapping wall | decision | bracket cost | controller CPU | comparison |
  | ---: | ---: | --- | --- | ---: | --- |
  | 8 | 13.620, 14.877 s; median 14.248 s | 4 -> 1 -> 4 -> 5; final 5 | 2 points, 16 samples, 1.353--1.453 s | 0.404--0.465 ms | median/oracle 1.0357 versus the 13.7575 s pin-5 median; <=5% gate passed |
  | 32 | 7.938 s | 4 -> 1; final 1 | 1 point, 2 samples, 0.800 s | 0.205 ms | faster than the reviewer's 8.61 s pin-1 median and inside the 1--3 gate |
  | 64 | 4.142, 4.618 s; median 4.380 s | 4 -> 1; final 1 | 1 point, 2 samples, 0.701--0.800 s | 0.216--0.308 ms | below the noisy same-binary pin-1 median of 4.925 s (3.231, 6.620 s) |

  At t8, median controller CPU was 0.434 ms, or 0.00305% of one core over
  mapping wall. Every adaptive and fixed output checked here retained canonical
  digest `293401b75cfd9a9003e3c83b0fe69070b68246f7c3b5a3ae69296674f36bb909`.
  Evidence is retained under
  `/scratch3/rob/tmp/tb-opening-bracket-final-t{8,32,64}-*`; the rejected 300 ms
  horizon trials are retained as t8 r3/r4 and showed why a favorable single run
  is not a sufficient gate.
- The remaining scATAC-specific numeric choices are deliberately separated by
  semantics: opening four is a leaveable performance hint, reader batch 1024 is
  drain/throughput geometry, progress publication 64 is measurement resolution,
  and the 5 s steady probe is the independently measured sparse-responsive
  policy. Only the opening and sparse policy enter the broker configuration;
  none is a safety floor. This does not literally reduce every scATAC numeric
  choice to two, because folding batching, measurement resolution, and settled
  responsiveness into one knob would recreate the category error in another
  form. Each surviving value has measured provenance in this ledger and the
  design document.
- Final verification passed `cargo fmt --all -- --check`, `git diff --check`,
  `cargo check --no-default-features`, strict all-target Clippy with and without
  `rapidgzip-cpu-accounting`, all 16 controller unit tests, all 44 synthetic
  control-law tests, all 243 runnable piscem library tests (one fixture test
  ignored), and all three local oracle-boundary tests. The final release binary
  was rebuilt after the last source change before the real-data runs.

### 2026-08-09 — Bracket-review bimodality correction

- `bracket-review.md` reproduced a real defect in `e78757c`: t32
  selected model producer one in 5/8 runs but exhausted the bracket and restored
  the unevidenced producer-four opening in 3/8. The 9.68 s median and
  7.97--12.77 s range showed why the earlier single t32 run was not a valid
  qualification result. The review's constant-count concern is explicitly
  closed with no action; its housekeeping warning was honored by keeping
  `.gitignore`, `scripts/`, manifests, data, and result trees out of the commit.
- The first literal correction—requiring non-overlapping intervals and a
  greater-than-5% loss before rejecting the model—was necessary but not
  sufficient. A fresh eight-run gate still produced one separated false
  regression and restored producer four. Worse, two t8 guards treated their
  model/opening result as inconclusive, retained producer one, and slowed to
  18.95--19.23 s. These rejected experiments remain at
  `/scratch3/rob/tmp/tb-bracket-inconclusive-fix-t32-r*` and
  `/scratch3/rob/tmp/tb-bracket-confirmed-t{8,32}-r*`.
- The completed design separates trigger evidence from retained evidence:

  1. the short model point remains the common fast path;
  2. an apparent regression is extended to the ordinary ratification sample
     count before exploration;
  3. persistent slack/source cap evidence at or below the model target can
     independently confirm it and suppress redundant exploration;
  4. otherwise a regressed or inconclusive comparison may test adjacent points,
     but each must beat the higher measured opening/model rate by more than 5%
     with non-overlapping intervals;
  5. the model remains the fallback if no point wins or the deadline expires.

  The opening is therefore a candidate-placement pivot and high-water evidence
  bar, never a fallback. This is allocation-derived and generic: t32's measured
  useful cap one confirms model one, while t8 retains headroom and tests five.
  No modality, budget threshold, or new tuning constant was introduced.
- Tests pin the statistical rejection rule, independent cap confirmation,
  model-versus-opening fallback outcome, asymmetric confirmation extension,
  deadline fallback, and the paired opposite-shape behavior. Real evidence also
  showed the model deadline fallback directly.
- Final validation used release binary SHA-256
  `78cebc9c499e93ee93e96abe520a786cb9efd18238cfb742e89ec41ed6b03317`
  on the 2,000,000-read / 528,410-mapped fixture:

  | cell | n | final/outcome histogram | mapping seconds | bracket wall | controller CPU |
  | --- | ---: | --- | --- | --- | --- |
  | scATAC t32 | 8 | producer 1: 8/8; `model_selected`: 7, `budget_exhausted`: 1 | 4.526--11.859; median 7.808 | median 0.800 s; normal 0.700--1.001 s; max 4.000111 s | 0.129--0.363 ms; median 0.173 ms |
  | scATAC t8 | 2 | producer 5 / `alternative_selected`: 2/2 | 14.412--14.839; median 14.625 | 1.052--1.153 s | 0.296--0.371 ms; median 0.334 ms |

  The t32 median is 0.907 times the reviewer's 8.61 s pin-one reference and
  passes its <=1.05 gate; no run ended outside producer 1--3. The one deadline
  report exceeded the nominal four seconds by only 111 microseconds of scheduler
  wakeup latency and retained producer one rather than restoring four. Every
  adaptive output and both contemporaneous t8 pin-five controls produced
  canonical digest
  `293401b75cfd9a9003e3c83b0fe69070b68246f7c3b5a3ae69296674f36bb909`.
  Final evidence is at
  `/scratch3/rob/tmp/tb-bracket-cap-gate-final-t{8,32}-r*`; fixed t8 controls are
  in the corresponding `t8-pin5` paths.
- Real discrete-decision gates now require repetition count, full range,
  final-allocation and outcome histograms, and median. A previously bimodal cell
  requires at least eight post-fix repetitions, zero out-of-region outcomes, and
  a per-run bracket-wall bound. A single run is explicitly diagnostic only.
- Final verification passed formatting and diff checks, feature-off compilation,
  strict all-target Clippy with and without `rapidgzip-cpu-accounting`, 19
  controller unit tests, 45 control-law tests, 243 runnable piscem library tests
  (one ignored fixture), and the three local oracle-boundary tests.
