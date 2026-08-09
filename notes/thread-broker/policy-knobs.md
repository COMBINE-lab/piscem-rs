# Modality policy constants, and a proposal for self-bracketing openings

For the thread-broker owner. Two parts: a finding about the hard-coded scATAC
constants (§1–§4), and a proposed longer-term mechanism to retire most of them
(§5–§7). §8 is the working agreement.

Context: found while reviewing `f36def1` on
`fix/thread-broker-scatac-t32-floor`, after the five-modality oracle matrix on
merged `main` (`ed59850`).

**The headline is not "the number is wrong".** The number is fine today. The
finding is that a *performance belief was written into a safety constraint*, that
the belief was fitted to the edge of a grid that did not bracket the optimum, and
that the fix repeated the same procedure. That is a process defect and it will
recur until the mechanism changes.

---

## 1. What is hard-coded

`src/cli/map_scatac.rs`, after `f36def1`:

```rust
fn scatac_broker_config(budget: usize) -> BrokerConfig {
    let mut config = broker_config_for_budget(budget);
    config.steady_probe_interval = Some(Duration::from_secs(5));   // (a)
    if budget <= 8 {
        config.nonlinear_probes = true;                            // (b)
        config.nonlinear_probe_samples = 12;                       // (c)
        config.nonlinear_blackout_samples = 12;                    // (d)
    } else {
        config.min_producer_slots = 2;                             // (e)  CONSTRAINT
        config.nonlinear_probes = false;                           // (f)
    }
    config
}

fn scatac_initial_decode_slots(budget: usize) -> usize {
    let measured_opening = if budget <= 8 { 4 } else { 2 };        // (g)  HINT
    measured_opening.min(budget.saturating_sub(1)).max(1)
}
```

Plus, outside the broker config: `SCATAC_READER_BATCH_SIZE = 1024` (against
`READER_BATCH_SIZE = 16384`) and the 64-vs-256 record progress publisher with
its per-processor cache-padded shards.

Eight workload-specific constants, and the `budget <= 8` threshold appears four
times.

---

## 2. The category error

`(e)` and `(g)` encode the same belief — "about 2 slots is right at large
budgets" — but in two different kinds of slot, and only one of them is
appropriate.

`min_producer_slots` is a **safety** constraint. Its own documentation says so:

> *"`rapidgzip` used to read this limit while choosing a backend, so anything
> below two committed a decoder to sequential decoding irreversibly … Measured on
> the shared-pool decoder: an aggregate limit of 1 leaves the path at
> `DenseMembers` and still reports a demand of 12, so admission is genuinely
> decoupled and **the floor can be a single slot**."*

So the safety-justified value is **1**. Setting it to 2 to express a performance
opinion makes the controller **structurally unable to represent the model's own
answer**. At `t32` the model computes `ideal = 1`; `target = ideal.clamp(2, hi)`
= 2. It did not weigh 1 and reject it — 1 is unreachable, and no amount of
evidence can make it reachable.

A floor is a constraint the controller may never violate. An opening is a guess
it is free to leave. The belief belongs in the opening.

---

## 3. How it happened, twice

The normal-tier validation grid is `pin_fractions = [0.125, 0.25, 0.333]`:

| budget | grid |
|---|---|
| t8 | 1, 2, 3 |
| **t32** | **4, 8, 11** |
| t64 | 8, 16, 21 |

At `t32`, **4 is the smallest point in that grid.** It won, and the code recorded
*"Normal t32/t64 fixed surfaces repeatedly put producer 4 in the safe oracle
region"* — then wrote 4 into the floor, where it blocked everything below 4.

My denser grid was `{2, 4, 8, 11, 16, 24}`. **2 is the smallest point in that
grid.** It won, and 2 was written into the floor. Same procedure, same failure
mode. My grid did not sample 1 or 3 at `t32` either, so the new constant has
exactly the epistemic status of the old one. `f36def1`'s message half-concedes
it: *"At t64, producer 2 is within 4.3% of producer 1"* — 1 is known to be
competitive, and the floor forbids it.

### What it costs today: nothing measurable

Measured at `t32` on the 2 M-read fixture, two repetitions each:

| decode / map | median | runs |
|---|---:|---|
| 1 / 31 | 8.61 s | 8.88, 8.35 |
| 2 / 30 | 9.82 s | 10.34, 9.30 |
| 3 / 29 | 8.98 s | 8.56, 9.40 |
| 4 / 28 | 10.42 s | *(3 reps, matrix)* |

1, 2 and 3 are indistinguishable at this sample size — run-to-run spread is
±0.6 s on a 9 s run — and all clearly beat 4. **The current constant costs
approximately nothing.** That is the point: the value is not the problem, the
procedure that produced it is, and the procedure has now demonstrably shipped a
wrong constant once.

---

## 4. A mechanical check that would have caught both

> **If a policy constant is the smallest or largest point in the grid used to
> justify it, the grid did not bracket the optimum and the constant is
> unproven.**

True of `4` in `{4, 8, 11}`. True of `2` in `{2, 4, 8, 11, 16, 24}`. It is cheap
to assert in the oracle runner — refuse to report an "oracle" when the winning
point is on the boundary of the swept set — and it catches this exact mistake
without anyone having to remember it.

**Immediate suggested change** (three lines, no behaviour change when 2 is
right): drop `config.min_producer_slots = 2`, leaving the crate default of 1, and
keep the opening at 2 in `scatac_initial_decode_slots`. The run still *starts* at
2, so there is no convergence cost in the common case, but the model can solve to
1 and ratification decides on measurement rather than on a constant.

---

## 5. Why this strains the crate's generality

The inventory is lopsided, and it is worth being precise, because the crate is
**not** generally knob-heavy:

| modality | broker knobs set |
|---|---|
| bulk SE | none — shared `broker_config_for_budget` + the 3 universal builder calls |
| bulk PE | none |
| scRNA | none |
| Flex | none |
| **scATAC** | **8**, budget-branched, plus a bespoke progress-shard path |

Four of five modalities need zero tuning. One needs eight constants. So the
problem is not that the interface demands tuning in general — it is that scATAC
has physics the model cannot see, and the policy layer is carrying it.

**Specifically: allocation-dependent consumer scaling.** At `t8`, adding mapping
threads makes mapping *slower* — the true peak is 5 producer slots against a
model answer of 1, and going to 1 costs 41%. The busy-ratio cost model assumes
each stage's per-item cost is independent of its allocation, and here that is
false. Every one of the eight constants is a human-measured stand-in for a
measurement the controller is not currently able to make.

That matters beyond scATAC. The second adopter is salmon, and the honest warning
is: *if your workload has allocation-dependent stage costs, the model will not
see it, and you will end up writing constants like these — each valid only for
the datasets you measured.*

---

## 6. Proposal: an opt-in self-bracketing opening

### 6.1 The cost position, stated explicitly

**A bounded startup cost is acceptable and expected.** The purpose of the broker
is near-optimal allocation on long runs; spending a little at initialisation to
avoid being wrong for the whole run is a good trade.

**What is not acceptable is making the broker substantially heavier overall.**
Two hard constraints on any design here:

1. **Steady state must not change.** No increase in the settled sampling cadence,
   no additional recurring work, no change to the ≤0.001-core administrative CPU
   result. The bracket is a startup phase that ends.
2. **The startup cost must be explicitly budgeted, bounded, and reported.** A
   fixed cap in wall time and sample count, defaulted conservatively, surfaced in
   telemetry so it can be attributed rather than inferred.

Under the rule that runs shorter than ~30–60 s do not gate, a bounded startup
cost of order a second is by construction a shrinking fraction of any run the
broker exists to optimise.

### 6.2 The key property: cost proportional to disagreement

The bracket should run **only when the model's answer differs from the opening**.
If they agree — the common case for bulk, scRNA and Flex, which already open near
their optimum — it costs exactly zero. Only a workload where the model and the
opening disagree pays, and that disagreement is precisely the signal that one of
them is wrong.

This is what keeps the mechanism from being a tax on the four modalities that do
not need it.

### 6.3 Shape

Extend the existing policy enum rather than adding a parallel mechanism:

```rust
pub enum OpeningPolicy {
    /// Today's behaviour: trust the caller's opening and the cost model.
    Fixed,
    /// Bracket the opening at startup, within a hard budget.
    Bracket {
        max_points: usize,        // e.g. 3
        horizon: Duration,        // per point, e.g. 300 ms
        total_budget: Duration,   // hard cap across the whole bracket
    },
}
```

Reuse the machinery that already exists — it is most of the work:

- bounded point sequence with a tried-set,
- ratification with interval comparison and explicit inconclusive outcomes,
- retained-best with physical restore before the next probe,
- interior refinement for a non-adjacent retained/failing pair,
- independent confirmation horizons.

Two changes are needed:

1. **It must be able to bracket downward.** `nonlinear_probe_targets` currently
   only grows from the opening (`4 → 6 → 9 …`). A bracket must test the model's
   answer and the neighbours *between* it and the opening. This is also the
   direct fix for the `t32` class of failure.
2. **It runs at startup, keyed on model/opening disagreement**, rather than only
   after the model reaches the producer floor.

The model still leads. The bracket is a *confirmation*, not a sweep: test the
model's answer, and if it loses, one or two points between it and the opening.
Typically ≤3 points, and the `t8` trajectory under `f36def1`
(`4 → 1 → 4 → 6 → 7 → 6 → 5`) shows the machinery already converges in that
neighbourhood — it just cannot start from the right question.

### 6.4 Telemetry

Report `opening_bracket { points_measured, wall_nanos, samples, outcome }` so the
cost is auditable and the oracle runner can attribute it, exactly as
`nonlinear_probes` and administrative CPU are today.

---

## 7. Acceptance criteria

A design is done when all of these hold:

1. **scATAC `t8`** still selects producer 5 (dense-grid peak 13.76 s), within 5%.
2. **scATAC `t32`** selects 1–3 (all within noise of 8.6–9.8 s), within 5% —
   *without* `min_producer_slots` being used to encode it.
3. **bulk SE/PE, scRNA, Flex** are unchanged: no bracket triggered, or a measured
   cost inside the stated budget. These already land within ±8%; none should
   regress.
4. **Steady state is untouched** — settled cadence, sample counts and
   administrative CPU match the current numbers.
5. **Startup cost is bounded and reported**, and stays inside its configured cap
   on every cell.
6. **The scATAC constant count drops from 8 to ≤2** — ideally an opening hint and
   the reader batch geometry, with the rest derived at runtime.
7. **Every surviving constant carries provenance** (grid, dataset, date, measured
   optimum) and passes the §4 boundary check.

A useful sanity property for the test suite: the `t8` and `t32` regression tests
must not be simultaneously satisfiable by any single hard-coded opening. They
have opposite answers from identical model output, which is what makes them a
real pair.

---

## 8. Working agreement

Please work on your own branch off `dev`, not on `dev` directly.

- `main` and `dev` are both at `ed59850` and in sync with origin.
- `fix/thread-broker-scatac-t32-floor` (`f36def1`) is local-only, one commit
  ahead, and reviewed — the `t8` excursion cost it introduces
  (13.71 s → 14.52 s, from the added floor validation) has been accepted.
- Merge back into `dev` once it passes; `dev` is where this work accumulates
  until we are ready to release.

Release remains blocked independently on the two unreleased git-pinned
dependencies.

---

## Appendix: evidence

| what | where |
|---|---|
| `t32` low-end sweep (1, 2, 3) | `/scratch3/rob/tmp/tb-lowend-p{1,2,3}-r{1,2}/` |
| `t32` + scRNA oracle fill | `/scratch3/rob/tmp/thread-broker-oracle-fill/` |
| full 7-point `t8` sweep | `/scratch3/rob/tmp/thread-broker-oracle-t8-fine/` |
| `t8`/`t64` oracle matrix | `/scratch3/rob/tmp/thread-broker-oracle-t8-t64-merged/` |
| post-fix verification runs | `/scratch3/rob/tmp/tb-fix-t{8,32}-r*/` |
| the original failure write-up | `scatac-t32-failure.md` (§4 carries a correction) |

Binaries: `5853003e…` at `ed59850` for the matrices, rebuilt at `f36def1` for the
verification runs. All runs `--features rapidgzip-cpu-accounting`, 2,000,000
reads / 528,410 mapped, canonical output equality throughout.
