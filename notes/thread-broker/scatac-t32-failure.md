# scATAC `-t 32`: the controller holds a split above its own measured cap

Investigation handoff. Found while running the five-modality oracle matrix on
merged `main` (`ed59850`), binary built with `--features rapidgzip-cpu-accounting`.

**Severity: the single real performance failure in a 320-run matrix.** 25% slower
than the best split, reproducible in 3/3 runs, deterministic. Everything else in
the matrix landed within ±8% of optimal.

**Short version.** The cost model is *right* — it reports decode at 0.09% of
per-read cost and asks for the producer floor. A guard added for the scATAC `t8`
case discards that answer, and the mechanism meant to replace it (the bounded
nonlinear probe) never runs. The controller parks at its arbitrary opening split
of 4 for the entire job, **above its own measured `useful_cap` of 2**, and never
moves.

---

## 1. Reproduction

~13 s. Deterministic; 3/3 matrix reps and 2/2 manual runs behaved identically.

```bash
cargo build --release --features rapidgzip-cpu-accounting

RUST_LOG=thread_broker=trace target/release/piscem-rs map-scatac \
  -i /scratch1/rob/atac_test/chr12_idx \
  -1 /scratch3/rob/tmp/tb-atac-10x-R1.fq.gz \
  -b /scratch3/rob/tmp/tb-atac-10x-R2.fq.gz \
  -2 /scratch3/rob/tmp/tb-atac-10x-R3.fq.gz \
  -o /tmp/out -t 32 --quiet --dict auto
```

Compare against the optimum:

```bash
… --decoder parallel=… # aggregate pin; see PISCEM_DECODE_SLOTS=2 for the peak
```

---

## 2. Ground truth: the surface

Three repetitions per point, medians, same binary, counterbalanced blocks.

| decode / map | mapping |
|---|---:|
| **2 / 30** | **9.22 s** ← true peak |
| 4 / 28 | 10.42 s |
| 8 / 24 | 10.38 s |
| 11 / 21 | 13.40 s |
| 16 / 16 | 14.77 s |
| 24 / 8 | 17.23 s |
| *adaptive* | *11.51 s* |

`auto / peak = 1.249`. CPU ratio 1.169.

Note the surface confirms the model: decode slots are nearly worthless here —
24 of them costs 87% over the peak. The right answer is at or near the floor.

---

## 3. Observed controller behaviour

Telemetry, all three matrix reps:

```
rep1: traj=[4,4]  moves=0 reverts=0 resurveys=0  nonlinear_probes=0  override=false
rep2: traj=[4]    moves=0 reverts=0 resurveys=0  nonlinear_probes=0  override=false
rep3: traj=[4]    moves=0 reverts=0 resurveys=0  nonlinear_probes=0  override=false

model: share=0.0007–0.0013  ideal_producer_slots=1  effective_deadband=1
       useful_cap=2  useful_cap_reason=slack
plan:  decode_slots=4 (opening)      final_phase=steady   terminal_error=none
       pressure_vetoed_shrinks=0
```

Trace (default sparse-responsive cadence):

```
02:00:38.443  thread execution plan  effective_budget=32 mapping_threads=28 decode_slots=4
02:00:38.843  survey  idle 15.1%  pressure Satisfied  consumer 28 (live 28) producer 4 (active 0)
02:00:38.943  survey  idle  0.0%  pressure Satisfied  consumer 28 (live 28) producer 4 (active 0)
02:00:39.044  steady  idle  0.0%  pressure Satisfied  consumer 28 (live 28) producer 4 (active 0)
…            steady … producer 4 (active 0)   [for the rest of the run]
02:00:51.534  thread broker: settled at 28 mapping / 4 decode (0 moves, converged true, 300ms)
02:00:51.534  thread broker: decode is 0% of per-read cost -> wanted 1 slots, usable ceiling 2
```

Two things to notice:

* **`producer 4 (active 0)` for the entire run.** The decode side never occupies
  a single execution slot. The four reserved slots do nothing at all.
* It reports `converged true` after 300 ms. By its own definition it has
  succeeded, while sitting two slots above its own cap.

---

## 4. Root cause chain

> **Correction (2026-08-09).** An earlier revision of this document blamed the
> nonlinear floor guard in `controller.rs`. That was wrong, and the error is
> worth recording because the telemetry invites it: the report field
> `ideal_producer_slots = 1` is the **pre-clamp** model answer, not the target
> the controller acted on. I read it as the target. The real cause is below and
> lives in the application, not the control law.

`src/cli/map_scatac.rs`, `scatac_broker_config`, at the reviewed commit
`ed59850`:

```rust
if budget <= 8 {
    …
    config.nonlinear_probes = true;
} else {
    // "Normal t32/t64 fixed surfaces repeatedly put producer 4 in the safe
    //  oracle region … Preserve four as the application-measured floor/opening"
    config.min_producer_slots = 4;
    config.nonlinear_probes = false;
}
```

At `t32` that makes the chain trivial:

1. `Costs::solve` computes `ideal = 1` from the 0.09% cost share — correct.
2. `let lo = cfg.min_producer_slots.max(1)` → **4**, so
   `target = ideal.clamp(lo, hi).min(useful.max(lo))` = **4**.
3. `distance = target.abs_diff(split.1)` = `|4 − 4|` = **0**, which is below the
   deadband, so the Survey arm settles immediately. Hence `moves = 0` and
   convergence reported at 300 ms.
4. `nonlinear_probes = false` for budgets > 8, so `nonlinear_probe_enabled` is
   false and `nonlinear_probe_complete` starts `true`. Hence
   `nonlinear_probes = 0` — the probe was **switched off**, not failing to run.

**The control law behaved exactly as configured.** The defect is the
application policy: scATAC pinned both its floor and its opening at 4 slots for
every budget above 8, on the strength of an earlier sweep whose grid did not
include 2. The denser grid used here (`2, 4, 8, 11, 16, 24`) puts the optimum at
2, which the old policy had made unreachable — the floor was *above* the answer.

This is the coarse-grid problem showing up as a shipped default rather than as a
misleading measurement: a policy constant derived from a grid that never sampled
the optimum.

## 5. Why the obvious fixes are wrong

This is the hard part, and it is why the guard exists. **scATAC at `t8` produces
an identical model and needs the opposite answer.**

| | `t8` | `t32` |
|---|---|---|
| model `share` | ~0.001 | ~0.001 |
| model `ideal` | 1 | 1 |
| `useful_cap` | 2 (`slack`) | 2 (`slack`) |
| opening split | 4 | 4 |
| **true peak** | **5** | **2** |
| nonlinear probes | 3 | **0** |
| trajectory | `4→6→7→5` | `4` |
| result | **0.996 of peak** | **1.249 of peak** |

At `t8` the model and the cap are both badly wrong — the true optimum is 5, more
than double the cap — and the nonlinear probe is the only thing that finds it.
At `t32` the model is essentially right and the probe is unnecessary.

So these are all wrong:

* **"Honour the model when it says floor."** Breaks `t8`, where the floor answer
  costs 41% (19.36 s at 1 vs 13.76 s at 5).
* **"Never hold a split above `useful_cap`."** My first instinct, and wrong for
  the same reason: at `t8` the cap is also 2 while the optimum is 5. The cap is
  derived from producer slack and is not trustworthy when the producer is idle
  because the *consumer* is the constraint.
* **"Gate the guard on budget."** A band-aid; the discriminator is not the
  budget, it is whether allocation-dependent consumer scaling actually exists on
  this workload.

The real discriminator between the two cases is that at `t8` shrinking the
producer *hurts* (consumer scaling is negative — more mapping threads make
mapping slower), and at `t32` it *helps*. Neither the cost share nor the cap can
see that. Only measuring can, which is precisely what the nonlinear probe is for.

**So the fix is almost certainly not to weaken the guard, but to guarantee the
probe runs whenever the guard fires** — and to make the guard time-bounded, so
that if the probe has not produced a proven point within some horizon, the
model's answer is honoured rather than the arbitrary opening split being kept
forever.

Note the opening split has no evidential status here. `decode_slots = 4` comes
from the generic `budget / 8` heuristic, not from measurement. The guard's
comment calls it "the caller's measured safe opening", which is true at `t8`
(where 4 was empirically chosen) but false at `t32`.

---

## 6. Secondary finding: same-split overhead

Adaptive at 4 decode slots runs **11.51 s**; pinned at the same 4 slots runs
**10.42 s** — 10.5% slower at an identical allocation. Some of the 1.249 is
therefore not the split choice at all.

Candidates: the adaptive-only 64-record progress publisher (pinned paths use
256), or the pool's immutable `workers` ceiling differing between adaptive and
pinned admission. Worth isolating separately; it may also explain the opposite
sign seen on Flex `t8`, where adaptive was ~3% *faster* than the pin at the same
split.

---

## 7. Suggested regression test

Deterministic, in `control_law.rs`, using the coupled fake:

* budget 32, opening producer 4, producer cost ~0.1% of total, consumer scaling
  linear (no negative scaling), true optimum at the producer floor;
* assert the controller reaches ≤2 producer slots within N samples;
* keep the existing `t8`-shaped test (negative consumer scaling, true optimum
  above the opening) passing unchanged.

The pair is the point: one test must not be satisfiable by breaking the other.

---

## 8. Evidence

| what | where |
|---|---|
| oracle matrix, `t8`/`t64` | `/scratch3/rob/tmp/thread-broker-oracle-t8-t64-merged/` |
| full 7-point `t8` sweep | `/scratch3/rob/tmp/thread-broker-oracle-t8-fine/` |
| `t32` + scRNA fill (this cell) | `/scratch3/rob/tmp/thread-broker-oracle-fill/` |
| sparse-cadence trace | `/scratch3/rob/tmp/tb-diag/trace.log` |
| dense-cadence run (`…PROBE_INTERVAL_MS=25`) | `/scratch3/rob/tmp/tb-d3/out/map_info.json` |
| manifests | `scripts/thread_broker_oracle.{t8-t64-merged,t8-fine,fill}.json` |

Binary `5853003ead3f0048…`, commit `ed59850`, `rapidgzip-core`
`276a41f7…`, `paraseq` `40ae74fe…`. All 320 matrix runs exited 0 with canonical
output equality.

---

## 9. Two smaller items from the same matrix

Not part of this failure, but found alongside it and cheap to fix:

1. **A no-op resurvey resets `time_to_converge`.** Flex `t8` rep 1 resurveyed
   once, made zero moves, and reported convergence at 17.24 s instead of
   0.075 s — failing the convergence gate on a run whose split never changed.
   Convergence should be judged on split stability, not on whether the decision
   was re-opened. This will keep producing false gate failures.
2. **`resurvey_distance: 1` is sensitive at small budgets** — that same spurious
   resurvey, 1 of 3 reps, on a stable workload.
