# Thread-broker review 2 — post-audit implementation

Review of `53388a6` (`feat(io): complete gated thread-broker implementation`) on
`feat/thread-broker`, against `audit.md`,
`completion-ledger.md`, and the code.

**Overall: the direction is right and the work is high quality.** One defect was
found that could destroy the output of a completed run; it is fixed in this
working tree (§1). Everything else recorded here is either validation that has
not been run yet or a design question worth an explicit answer — none of it
blocks merging the branch, and most of it is already tracked in the completion
ledger.

The ledger's own discipline deserves saying plainly: failures are retained with
their numbers, a mis-counterbalanced experiment was invalidated and rerun, and a
previously-claimed pass was withdrawn. That is the main reason the unverified
parts of it are credible.

---

## 1. Fixed in this tree — broker failure destroyed the output of a successful run

**Severity: high. Reachable. Silent.**

### What it was

All three CLIs propagated a broker failure out of the mapping function:

```rust
// map_bulk.rs:541, map_scrna.rs:621, map_scatac.rs:603
let mut r = broker
    .finish()
    .map_err(|e| anyhow::anyhow!("thread broker failed: {e}"))?;
...
result?;   // the actual mapping result was only checked after this
```

`run_bulk_pipeline` is invoked with `?` at `map_bulk.rs:263`. The RAD
chunk-count backpatch is at `:288-293` and `write_map_info` at `:317`. So on a
broker error, after mapping had **already completed successfully**:

1. `?` unwound past the backpatch;
2. the `num_chunks` placeholder written at `rad.rs:208` is **`0`**;
3. `map_info.json` was never written.

The run left behind a RAD file that opens cleanly and declares **zero chunks** —
an empty result from a job that succeeded, discovered only downstream.

### Why it mattered

None of the four `BrokerErrorKind` variants — `ResizeRefused`,
`ResizeTimedOut`, `ThreadSpawn`, `ThreadPanicked` — has any bearing on output
correctness. They are failures of *thread allocation*. The reads that were
mapped were mapped correctly and the split in force when the broker stopped
stays in force, under budget, for the rest of the run.

And `ResizeTimedOut` is reachable: it fires when a shrink is not acknowledged
within `resize_timeout`, and the ledger records that the scATAC reader batch was
cut from 16384 to 1024 records *specifically* because acknowledgement was not
completing inside the 2 s window. A slower disk, a larger batch, or a loaded
machine reaches it again.

Secondary effect: because `result?` came afterwards, a broker error also masked a
genuine mapping error — an I/O fault would have been reported as "thread broker
failed".

### The fix applied

New `crate::io::broker::settle(RunningBroker) -> Option<BrokerReport>`, with the
rationale written on it. It logs the failure at `warn!`, keeps the partial report
the error already carries — that report still describes what the broker did
before it stopped, which is what a reader of `map_info.json` needs — and returns
control to the mapping result.

Call sites became:

```rust
let broker_report = broker.and_then(crate::io::broker::settle).map(|mut r| { ... r });
```

The same principle was applied to `start()`: a refused initial resize or a
sampler that will not spawn now warns and maps **without** adaptation, at the
plan's startup split, which is already valid and under budget. Genuine
configuration errors (`build()`, `PISCEM_THREAD_BROKER_*` parsing) still fail
fast, correctly — those are user input, reported before any work is done.

No new test was needed for the invariant `settle` depends on: `control_law.rs`
already asserts `error.report()` is populated for runtime failures (`:1087`,
`:1190`).

### What to take from it, beyond the fix

The audit's TB-03 ("surface actuation failures") was right, and the
implementation followed it — but "surface" was implemented as "propagate". For
an advisory subsystem those are different things. **Please re-scan for any other
path where a tuning-subsystem failure can abort or corrupt work that already
succeeded.** The specific shape to look for is a `?` on broker/measurement
machinery that sits between the end of mapping and output finalisation.

---

## 2. Not started — the robustness matrix (`Milestone 4`, item 4)

Cold-cache, **FIFO**, mixed-input, source-bound, and inelastic coverage is the
one item still marked `[ ]`.

FIFO is the one to prioritise, and it is a correctness concern rather than a
performance one. The non-rewindable-input guard exists to prevent *data loss*: a
FIFO consumed by a probe cannot be re-read, and `open_gz_rapidgzip`'s magic sniff
plus re-open will hang forever on a pipe whose writer has exited (this is why
`-r <(zcat ...)` hung).

The 14 unit tests in `src/io/calibrate.rs` (including `non_rewindable_input_tests`
with its writer-less FIFOs) survive and pass, so the guard's logic is intact. But
**`ExecutionPlan` rewrote budget planning underneath it**, including
`reconcile_parallel_decoders`, which now retroactively converts a run to
`Serial` when no decoder handles open. That interaction has never been exercised
end to end.

Requested: an end-to-end cell for each of
- FIFO-only inputs,
- a mixed set (regular pair + FIFO pair) — the regular files must still get the
  parallel decoder and the FIFO must not be read before the run,
- a split group (regular R1, FIFO R2),

each asserting byte-identical output against the same reads as regular files, no
FIFO read before the run, and the right message at the right level for `auto`
versus explicit `parallel`.

---

## 3. Open validation, as the ledger already records

Listed here only so the handoff is complete; the ledger tracks all of these
honestly and none is hidden inside a completed milestone.

- **Gate H** — the five-modality × `-t {8,32,64}` matrix. Currently two
  modalities at a few budgets.
- **Gate E** — the real input-archetype matrix for producer measurement.
- **Gate F / TB-06** — five pinned response points per real cell; the five-point
  data that exists already shows completed-worker CPU spanning 288–398 ms for
  bulk and 320–404 ms for PBMC, which is material allocation dependence.
- **Gate G / TB-07** — the alternating-compressibility and lull fixture.
- **TB-08** — 1,000-trace ratification replay using variance measured from real
  runs rather than synthetic noise.
- **TB-12** — the `{0.25,0.5,1,2,5}` s short-run matrix and the 30-block
  duration matrix.

### 3.1 The unresolved overhead result deserves a decision, not just more data

The whole-broker gate failed on the 3.8 s bulk cell (mapping-wall upper bounds
1.663%, 1.558%, 1.744%), and the ledger correctly identifies the cause: a fixed
**36–40 ms convergence penalty**, not controller CPU — administrative CPU came in
at 0.0003 core, three orders of magnitude inside its gate.

A fixed startup cost is the *expected* shape for this design and it dilutes with
run length. The open question is not whether it can be measured away but whether
the universal ≤1% bound is the right gate for a term that is constant in absolute
time. Consider stating the criterion as an absolute convergence budget (e.g.
"≤50 ms fixed, and ≤1% on runs longer than N seconds") rather than a single ratio
that any sufficiently short job must fail. Either answer is defensible; leaving it
as a standing failure against an unreachable bar is not.

---

## 4. Design questions worth an explicit answer

### 4.1 The generic/specific boundary, and adopter number two

The crate is now **3,519 lines with 36 config fields and 37 builder methods**. I
checked: there is no modality logic inside it, and the nonlinear probe is
correctly factored — off by default, documented as an application policy, enabled
at the scATAC call site, skipped by freeze. That is the right structure.

But the stated purpose is shared infrastructure for salmon, cuttlefish 3 and
ruSTAR, and making scATAC work needed:

- four tuned knobs in `scatac_broker_config`;
- a bespoke cache-padded per-processor progress shard in `processors.rs`;
- a separate publication cadence for progress (64) from busy time (256);
- a modality-specific reader batch size (`SCATAC_READER_BATCH_SIZE`).

**Request: document the minimum viable adoption path.** Which of the 36 fields
must a new consumer think about, and which are defaults it can ignore? If the
honest answer is that a new adopter must rediscover the progress-granularity and
batch-geometry interactions from scratch, that is worth knowing *before* salmon
starts, and may argue for pulling some of it into the crate as a documented
`Consumer` implementation guide or helper.

### 4.2 The `Starved` shrink veto is an exception to a stated invariant

`ProducerPressure`'s documentation says pressure "may veto growth but must never
size the split", and that restriction carries the whole §2 argument about
one-sidedness. The new bounded correction at `controller.rs:1154-1166` lets a
`Starved` producer veto a model-requested *shrink*.

I think this is defensible — it is one-directional, it cannot choose a larger
target, it is counted in `pressure_vetoed_shrinks`, and it has a regression test
with 10× producer under-reporting. But the invariant text as written no longer
describes the code. Please update the doc comment on `ProducerPressure` to state
the exception and why it does not reopen the failure mode, so the next reader
does not have to reconstruct it.

### 4.3 The scATAC reader batch change reaches beyond adaptive runs

`SCATAC_READER_BATCH_SIZE = 1024` (from 16384) applies to scATAC regardless of
whether the broker is running. It was introduced to make consumer-shrink
acknowledgement complete inside the resize timeout — a broker concern — and gated
with a single fixed six-slot comparison (3.116 s → 2.851 s).

Two things to confirm: that serial and `--decoder parallel=N` scATAC runs do not
regress in wall time or peak RSS at this batch size, and that 16× more batches
does not change memory behaviour on large inputs. A 16× reduction in batch size
is a large change to justify from one comparison.

Related, and already noted as a follow-up in the ledger: pinned/no-broker scATAC
should publish progress at 256 rather than 64, since no controller exists to
switch phases.

### 4.4 `cargo clippy --all-targets` does not pass in a tree with the local assets

`examples/canonical_rad.rs` (gitignored, local-only) fails to compile under the
default feature set — it imports `read_*_rad_records`, which are configured out
without `parity-test`.

The ledger's pre-commit entry claims Clippy for "every thread-broker target and
for piscem library and binaries", which is accurate and narrower than
`--all-targets`. Both of those do pass; I verified them. But the standard command
a reviewer reaches for fails in any tree that has the validation harness, which
is a papercut worth removing: feature-gate the example, or document the required
feature next to it.

---

## 5. What I verified independently

| check | result |
|---|---|
| `cargo fmt --all -- --check` | passes |
| `cargo clippy -p thread-broker --all-targets -- -D warnings` | clean |
| `cargo clippy --lib --bins --features rapidgzip -- -D warnings` | clean |
| `cargo clippy --lib --bins -- -D warnings` | clean |
| `cargo test -p thread-broker --all-targets` | 10 unit + 39 integration pass |
| `cargo test -p piscem-rs --features rapidgzip --lib` | 236 pass, 1 pre-existing ignored |
| `cargo check --features rapidgzip` after the §1 fix | clean |

Only the pre-existing `proc-macro-error2 v2.0.1` future-incompatibility note from
a third-party dependency remains, in every configuration.

## 6. What I did not verify

Everything measured. All timing, overhead, oracle-gap, convergence, CPU, RSS and
canonical-digest numbers in the ledger are taken as reported; I re-ran none of
them, and the raw evidence lives in ignored `/tmp` paths that I did not inspect.
The §1 defect was found by reading, not by reproducing it.

---

## 7. Two things the audit got right that are worth preserving

Recording these so they are not lost in a list of open items.

**TB-06 corrected a real error in the original design.** The design asserted that
the cost share *should* be budget-invariant and treated the measured 37% at
`-t 32` versus 26% at `-t 64` as evidence of an instrumentation bug. That was too
strong: allocation-dependent service cost is physical — cache and
memory-bandwidth contention — so the discrepancy may be entirely real. Replacing
cross-budget invariance with same-split repeatability plus agreement against
cumulative CPU ground truth is a strictly better instrumentation gate.

**The scATAC failure was found by running the case the design flagged as
untested.** Auto collapsing to one decode slot — 16 s against a 3 s pinned oracle
— is exactly the linearity assumption of the closed form breaking under negative
mapping-thread scaling. Finding it required running a modality nobody had run,
which is precisely the gap the original document listed under "where I think this
is weakest". The response (bounded geometric probe, then interior refinement,
with retained-best restore and independent confirmation horizons) is the right
shape and is well tested.

---

# Addendum — what "TB08 needs per-modality epochs" means, and how to run it

Added after the response in `review-2-response.md`. The response's
handling of TB08 is correct and I am not asking for it to be redone; this is the
concrete gate design that the phrase was standing in for.

## What TB08 is actually testing

The `Ratify` phase decides keep-or-revert by comparing pipeline throughput before
a move against throughput after it. That is a statistical hypothesis test on
noisy samples, and TB08 asks whether the test works: does it catch a real
regression, and does it avoid reverting a good move on noise?

## What the real replay found, and why it is not interpretable

- synthetic, 1,000 seeds at 2% window noise: ≥95% detection, <1% false confident
  rejection — passes;
- pooled real traces: 10 windows detected a true 10% regression **77.4%** of the
  time; 20 windows detected **3.8%**, with the population CV rising to 38.9%.

**Doubling the sample size made detection roughly twenty times worse.** That
cannot happen when sampling one population — more data always helps. It is the
signature of pooling several populations. The pooled variance is then dominated
by *between-regime* differences rather than within-regime noise, so each extra
window adds regime variance faster than `√n` shrinks the standard error.

So the experiment as run asks: "can you detect a 10% shift in a mixture of
startup ramp, post-move transient, steady state, burst phases and several
modalities?" The controller never faces that question. It compares one bounded
post-blackout window at split A against another at split B, inside a single
stationary stretch. A failure of the pooled experiment tells us the mixture is
heterogeneous — which we already knew — not whether ratification works.

## What an epoch is here

A **stationary span**: a contiguous stretch of one run, one modality, one budget,
one split, after blackout, during which the workload's cost structure is not
changing. Within-epoch window-to-window variance is the quantity ratification
actually contends with.

**The controller already emits the boundaries.** `epoch` increments only on a
detected persistent workload change (`controller.rs:1232`), is cleared alongside
the rejection set, and is already serialized as `final_epoch` and on each
rejection record. The strata are labelled; the gate simply has to stop pooling
across them.

## Gate design

1. Emit or parse per-window throughput tagged with
   `(modality, budget, epoch, phase, split)`.
2. **Discard windows in warm-up, blackout and draining.** Those are transients by
   construction and the controller never ratifies on them; including them is most
   of the pooled variance.
3. Estimate within-stratum CV per `(modality, epoch, split)` and report the
   distribution. If the design is sound this should sit far below the pooled
   38.9%; if some stratum is genuinely that noisy, that is a real and specific
   finding.
4. Replay 1,000 seeds **per stratum**, injecting a 10% regression, using that
   stratum's own measured variance.
5. Gate per stratum at the production default of 10 ratification windows:
   ≥95% detection, <1% false confident rejection on equal rates. **Report the
   worst stratum; do not average across them** — averaging is the same mistake
   as pooling, one level up.
6. If a stratum cannot reach 95% at 10 windows, the answer is a modality-specific
   `ratify_samples` for that stratum, not a global increase.
7. Keep the pooled result published as the information boundary it is.

## What not to do

Do not respond to the pooled number by sampling faster or publishing progress
more often. The response already says this and it is right; it is worth
underlining because the 3.8% figure looks alarming and the instinct is to fix it
with more data — which is exactly what produced it.
