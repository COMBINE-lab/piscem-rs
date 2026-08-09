# Response to thread-broker review 2

Review target: `53388a6` (`feat(io): complete gated thread-broker implementation`)

Remediation commit: `2b6ed60` (`feat: harden adaptive thread broker`)

Upstream rapidgzip dependency: `276a41f77fb927e24cb0898a638a08b21eb048c6`

Thank you for the careful review. We agree with the high-severity output-safety
finding and with the distinction between a merge-ready controller and complete
release qualification. The implementation findings have been addressed. The
remaining measurements are named below and have not been converted into implied
passes.

## Summary disposition

| review item | disposition |
| --- | --- |
| Broker failure could suppress RAD finalization | Fixed and covered by injected tests |
| FIFO/mixed/split-group robustness matrix | All requested end-to-end cells passed |
| Gate E producer measurement | Complete with exact upstream accounting and direct no-overhead gate |
| Gate F / allocation dependence | Bounded corrections implemented; selected normal cells pass; comprehensive response surfaces remain opt-in |
| Gate G flow/cap transients | Both persistent-source and temporary-lull matrices passed |
| TB08 ratification confidence | Synthetic gate passes; pooled bursty-real criterion is retained as an information-limit failure |
| TB12 overhead | Dual absolute/direct/fractional criteria implemented and policy costs measured; normal duration series remains unrun |
| Minimum adopter contract | Documented in crate-level rustdoc |
| `Starved` invariant mismatch | Documentation corrected and behavior regression-tested |
| scATAC batch/progress scope | Progress is policy-scoped; retaining the small batch is supported by serial and pinned measurements |
| Local all-target Clippy failure | Fixed; strict all-target Clippy passes in both feature configurations |

## 1. Broker failures are now advisory and cannot destroy completed output

We agree with the review's diagnosis and severity. The initial `settle` fix was
preserved and generalized into `AdvisoryBroker`, `BrokerDiagnostics`, and a
separate serializable `BrokerFailure` with stages for producer-measurement
startup, controller startup, and controller runtime.

The production behavior is now:

- user/configuration errors are rejected before mapping starts;
- measurement or controller startup failure logs a warning and leaves the safe
  planned split in force;
- resize refusal, timeout, counter regression, or controller panic stops further
  adaptation but does not determine the mapping result;
- any partial `BrokerReport` and the independent failure stage/message are
  retained in `map_info.json`;
- RAD finalization and `map_info.json` writing remain controlled by the mapping
  result, so a tuning failure cannot replace a real mapping error or suppress
  successful output.

We re-scanned the post-mapping path for the same failure shape. The broker finish
path no longer returns a propagated error. Fixed-measurement finishing is also a
diagnostic snapshot rather than a fallible post-mapping operation. Explicitly
requested validation instrumentation and invalid environment configuration
remain fatal before results are accepted, which is intentional.

Injected tests cover producer-measurement startup failure, initial resize
refusal, and controller panic. They assert that the failure is retained without
propagation. `io::map_info::tests::writes_broker_and_measurement_diagnostics`
checks the machine-readable failure serialization.

## 2. FIFO and mixed-input coverage is complete

The local one-shot harness ran the requested matrix under both `auto` and
explicit `parallel`:

- FIFO-only input;
- regular plus FIFO groups;
- regular R1 plus FIFO R2 in one split group.

All six treatments completed with no lingering or failed writer, named every
stream, used INFO for the automatic downgrade and WARN for an overridden
explicit request, matched regular-file counts, and produced byte-identical RAD
output. FIFO-only reconciled to serial `8/0`. Mixed and split-group runs retained
adaptive `6/2`, demonstrating that the regular gzip inputs still use the
parallel decoder without pre-reading the FIFO.

The harness is `scripts/thread_broker_fifo_matrix.py`; it and its generated
evidence are intentionally local and ignored. This suite is now a lifecycle and
input-path risk trigger rather than part of every routine controller iteration.
Cold-page-cache validation remains separate because we do not yet have a
portable, non-privileged way to guarantee a genuinely cold page cache.

## 3. Validation disposition

### Gate E — complete

The sampled producer-busy reconstruction was rejected after it showed 17--313%
lag error on sparse, bursty, and stored decoder paths. Production now reads
rapidgzip's feature-gated exact cumulative executing-region counter. No recurring
piscem sampler is started for this signal.

The upstream implementation uses a generation/transition snapshot protocol to
reject partial and ABA snapshots. With the feature disabled, the hot-path clock
reads, counter updates, transition operations, and related branches compile out.
Thirty counterbalanced direct-decoder pairs passed the no-measurable-overhead
gate with identical output:

- CPU ratio median `0.99577`, one-sided upper-95 `0.99858`;
- wall ratio median `0.99823`, one-sided upper-95 `1.00000`.

The final real Gate E matrix passed all 18 cells across six input archetypes and
three repetitions at
`/scratch3/rob/tmp/thread-broker-gate-e-real-epoch-v5-final`.

### Gate G — complete

The final persistent-source matrix passed all six cadence/repetition cells at
`/scratch3/rob/tmp/thread-broker-gate-g-real-v7-final`. Caps were identified in
`0.376--0.400 s`, exact paired producer-share bias was
`0.0003--0.0029` percentage points, occupancy stayed within budget, and every
canonical digest matched.

The separate temporary-lull recovery matrix passed all six cells at
`/scratch3/rob/tmp/thread-broker-gate-g-flex-recovery-v11-final`. Every eligible
event recovered its required allocation and cleared harmful cap evidence with
zero observed delay; paired share bias was `0.00024--0.0144` percentage points,
occupancy stayed within budget, and canonical output matched.

### Gate F / Gate H — controller evidence sufficient for landing; catalog not exhausted

The runner now enforces canonical content, budget and occupancy trajectories,
balanced mode positions, CPU, RSS, convergence, and an oracle comparison. The
final selected normal cells include:

- scATAC t32: adaptive median `13.594 s`, pin-4 oracle-region median
  `13.776 s`, adaptive/oracle upper-95 `1.0418`, CPU ratio `0.9868`, RSS ratio
  `0.9995`, canonical equality, and convergence;
- scRNA t32 after the cap-confirmation fix: adaptive `20.555 s`, pin-11
  `19.902 s`, oracle upper-95 `1.0760`, median one-third-baseline ratio
  `1.0328`, CPU ratio `0.9770`, RSS ratio `1.0156`, canonical equality, zero
  resurveys, and convergence in `0.3--2.1 s`.

Allocation dependence is treated as physical rather than automatically blamed
on instrumentation. The controller has bounded application-selected nonlinear
probing, local interior refinement, a `Starved` shrink veto, empirically
qualified modality openings, and conservative cap confidence. It does not fit
an unconstrained response curve online.

We are not claiming that the full five-modality/budget catalog ran on the final
binary. Completing the remaining normal cells is release qualification. The
540-run five-pin comprehensive matrix is opt-in for nightly runs or changes that
touch the controller/model assumptions.

### TB08 — the failure is retained, but it is not evidence to sample faster

The 1,000-seed synthetic ratification gate passes with explicit inconclusive
outcomes and isolated-zero protection. Replaying pooled final real traces does
not support the proposed 95% detection-power requirement: 10 samples detected a
true 10% regression in 77.4% of replays, while a 20-sample pooled run detected
only 3.8% as the heterogeneous population's CV rose to 38.9%. False confident
equal-rate rejection and isolated-zero rejection remained zero.

This is an applicability failure in the experiment, not a reason to increase
hot-path publication or controller frequency. Pooling startup, post-move,
modality, and burst regimes violates the stationary population assumed by the
test; simply adding windows can add regime variance faster than information.
The production default remains 10 ratification windows. A future real TB08 gate
must define stable per-modality epochs first and must preserve the failed pooled
result as an information boundary.

### TB12 — overhead criterion decided; long duration series still optional

We adopted the review's proposed distinction between fixed cost and fractional
cost. The validation runner reports mapping wall, process wall, aggregate CPU,
and direct controller/sampler CPU. Short startup-dominated cells use an absolute
`<=50 ms` statistical bound plus the direct administrative-CPU gate rather than
pretending a fixed 36--40 ms term can meet a universal fractional threshold.
The one-minute cell retains the fractional `<=1%` regression backstop as well as
the absolute bound. A normal-tier 10-minute cell is descriptive; 30-block
long-duration confidence belongs to the comprehensive tier.

The `1%` process-CPU criterion is only a broad backstop. At 64 busy threads it
could hide 0.64 core, so it is not the definition of cheap administration. The
primary direct gate is `<=0.001` core for responsive policies and `<=5 ms` total
controller/sampler CPU for freeze policies.

Measured policy costs on the 2-million-record, approximately 14--20 second
scATAC workload were:

| policy | mapping result versus pin-5 oracle | direct administrative CPU | interpretation |
| --- | ---: | ---: | --- |
| pinned/no broker | `1.000` | none | baseline |
| responsive, 25 ms | `1.217` | `4.965 ms`, `0.02936%` core | older search selected 6/2; allocation error dominates administration |
| sparse-responsive, 5 s | `1.208` | `1.408 ms`, `0.00835%` core | same split as responsive, 72% less administrative CPU |
| model-only freeze | `1.422` | `0.401 ms`, `0.00197%` core | near-free but wrong 1/7 split on this non-monotone surface |
| full-calibration freeze | `1.016` median / `1.037` upper-95 | `<=2.582 ms` | calibrated into the oracle region, then stopped recurring work |

A later eight-pair same-controller comparison isolated steady monitoring:
sparse/25-ms mapping was `0.9921/1.0210` median/upper-95, aggregate CPU was
`0.9941/1.0196`, and median administrative CPU fell from `4.902` to `2.623 ms`.
Sparse-responsive is therefore the closest policy to the original low-overhead
responsive design. Full-calibration freeze is the lower-recurring-cost choice
when the workload is stable and startup calibration can be amortized.

The 84-run normal 5-second/1-minute/descriptive-10-minute policy series was not
run before landing because it would add roughly 80--90 minutes and would not
exercise a new correctness path. It remains an explicit pre-release or
performance-change suite, not a hidden pass.

## 4. Design questions

### 4.1 Minimum viable adoption path

Crate-level rustdoc now has a `Minimum viable adoption` section. A new adopter
must provide only:

1. one truthful shared execution-slot budget, a valid opening, and real floors;
2. monotonic, non-blocking cumulative busy/progress counters;
3. truthful resize acknowledgement, including admitted work until retirement;
4. directional producer pressure rather than a desired thread count;
5. an explicit advisory-versus-fatal lifecycle policy.

Warm-up, smoothing, blackout, ratification, cap history, deadband, resurvey,
regression tolerance, and nonlinear probing stay at defaults until a named
measurement failure justifies changing them. The guide also names the second-
adopter qualification set: canonical equality, pinned response sweep,
source-bound and inelastic paths, a regime change, and both direct and whole-
process overhead measurements. Progress cadence and batch geometry remain
application concerns until repeated adopters demonstrate a genuinely generic
abstraction.

### 4.2 `Starved` shrink-veto invariant

The `Producer::pressure`, `ProducerPressure`, and `Starved` documentation now
states the exception. A persistent `Starved` signal may retain an occupied
current producer allocation when the cost model requests a shrink. It cannot
choose a target, request growth, or raise the current allocation. Source-bound
and inelastic evidence can cap growth; the cost-share model remains the only
sizing mechanism. The veto is reported and has a regression test with 10x
producer-cost under-reporting.

### 4.3 scATAC batch size and progress publication

Pipeline geometry is chosen after decoder-handle reconciliation and serialized
as `pipeline_tuning` in `map_info.json`. Fine 64-record completed-progress
publication is adaptive and decision-scoped; pinned, per-file-pinned, serial,
and reconciled-to-serial paths use the generic 256-record cadence. Busy-time
clock reads remain on their independent coarse cadence.

We measured rather than assumed whether the 1,024-record reader batch should be
adaptive-only. On the same two-million-record workload, 1K versus 16K batch
medians were:

| policy | mapping wall ratio | process CPU ratio | peak RSS ratio |
| --- | ---: | ---: | ---: |
| serial | `0.9839` | `0.9954` | `0.9389` |
| pinned | `0.9586` | `0.9690` | `0.9874` |

The smaller batch therefore had an independent performance/RSS benefit and
remains the scATAC default across policies. The result covers the representative
long fixture, but an even larger-memory scATAC rerun remains an appropriate
scATAC-hot-path risk trigger rather than a universal merge gate.

### 4.4 Standard all-target linting

The ignored canonical RAD helper now has a default-feature stub and gates its
parity-only imports. The final tree passes both:

```text
cargo clippy --all-targets -- -D warnings
cargo clippy --all-targets --features rapidgzip-cpu-accounting -- -D warnings
```

The only remaining message is Cargo's existing future-incompatibility notice
for third-party `proc-macro-error2 v2.0.1`.

## 5. Final verification and remaining catalog

The landing checkpoint passed:

- `cargo fmt --all -- --check`;
- `git diff --check`;
- strict all-target Clippy with and without rapidgzip CPU accounting;
- 15 controller unit tests and 41 control-law integration tests;
- 243 piscem library tests, with one fixture-dependent test ignored;
- feature-off compilation and the CPU-accounting feature configuration.

The remaining catalog is deliberately tiered:

- **light:** routine controller/unit feedback plus the common bulk-SE/Flex and
  short policy cells;
- **normal:** remaining five-modality `t8/t32/t64` oracle cells and the 84-run
  policy-duration series; this is the intended pre-release qualification;
- **comprehensive:** 540-run oracle and 360-run policy matrices, complete short-
  duration grid, cold-cache once portable, and risk-triggered Gate E/G,
  FIFO/mixed-input, scATAC geometry/publication, and rapidgzip direct A/B reruns.

Two non-controller release items remain: the git-pinned upstream dependencies
must be released or receive an explicit release exception. Rapidgzip's required
changes are merged upstream; the dependency is still unreleased.

## Conclusion

The review's implementation findings are resolved in `2b6ed60`. The controller
is suitable to land without waiting for the exhaustive matrices. That statement
does not claim complete release qualification: the normal duration series and
remaining modality cells are still named work, TB08's pooled-real failure is
preserved, and comprehensive reruns remain opt-in by risk.
