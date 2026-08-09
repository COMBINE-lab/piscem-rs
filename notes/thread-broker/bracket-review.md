# Review of `e78757c` — self-bracketing openings

Reviewing `e78757c` (`feat: self-bracket thread broker openings`) on `dev`,
against `policy-knobs.md` and
`policy-knobs-response.md`.

## Verdict

**Not ready — one specific defect. Everything else verifies.**

The design is right and the category error is properly fixed. But scATAC `t32`
is **bimodal**: 3 of 8 runs land on producer 4, the exact pre-fix wrong answer,
after exhausting the bracket's full time budget. As it stands the principled
mechanism is slightly *worse in median and much worse in variance* than the
hard-coded floor it replaces, on the cell that motivated it.

It is one small fix and one re-gate away.

---

## 1. What verifies

I reproduced the claimed validation exactly. The release binary built here has
SHA-256 `481b2be0163aa15c9feed1e1ac4f8c31d9f0c309447cf70809e3fa78f7380126` —
**bit-for-bit identical** to the one in the response — so these are measurements
of the same artifact, not a lookalike.

| claim | result |
|---|---|
| `cargo fmt --all -- --check` | pass |
| `git diff --check` | pass |
| `cargo check --no-default-features` | 0 errors |
| strict all-target Clippy, default | 0 errors |
| strict all-target Clippy, `rapidgzip-cpu-accounting` | 0 errors |
| 16 controller unit tests | 16 passed |
| 44 control-law tests | 44 passed |
| 243 piscem lib tests, 1 ignored | 243 passed, 1 ignored |

And the substantive claims:

* **Category error fixed.** The performance-derived `min_producer_slots` override
  is gone; the safety floor is 1 at every budget. `scatac_initial_decode_slots`
  is a uniform 4 with its provenance recorded in the comment, explicitly a
  leaveable hint.
* **Knob reduction is real.** `scatac_broker_config` is now three lines of
  configuration — sparse steady probe, opt-in bracket. The four nonlinear knobs
  and the `budget <= 8` policy branch are gone.
* **Shared config untouched.** The only change to `broker_config_for_budget` is
  two added *test assertions* that other modalities stay `OpeningPolicy::Fixed`.
  So the argument that bulk/scRNA/Flex need no new matrix is sound at the config
  level.
* **Cost proportional to disagreement — confirmed empirically.** Flex `t8`
  reports `outcome: not_configured, points_measured: 0, wall_nanos: 0`. Modalities
  that do not need the bracket pay literally nothing.
* **The controller refactor did not regress a `Fixed` modality.** This needed
  checking independently, because `controller.rs` changed 736 lines and the
  response's argument only covers configuration. Flex `t8`: **26.00 s / 26.94 s**
  at split 2 with 0 moves, against 25.77 s before the change — within noise, same
  split, bracket inert.
* **scATAC `t8`** — 2/2 runs select producer 5, 14.19 s / 14.44 s, matching the
  reported 14.25 s median.
* **scATAC `t64`** — **6/6** runs select producer 1, 4.08–5.48 s, all
  `model_selected` with 0.70–0.90 s bracket cost. Stable.

The bracket telemetry deserves specific credit: `outcome`, `points_measured`,
`samples` and `wall_nanos` are what made the defect below diagnosable in minutes
rather than requiring a debug build.

---

## 2. The blocker: scATAC `t32` is bimodal

The response reports `t32` from **a single run**, which caught the good mode.
Eight repetitions on the same binary and fixture:

| outcome | n | final split | mapping seconds |
|---|---:|---:|---|
| `model_selected` | 5 | **1** ✅ | 7.97, 8.12, 9.03, 9.27, 10.10 |
| `budget_exhausted` | 3 | **4** ❌ | 10.80, 11.47, 12.77 |

* median **9.68 s**, range **7.97–12.77 s** (60% spread)
* **37.5% of runs end on producer 4** — the pre-fix wrong answer
* every run produced 2,000,000 reads / 528,410 mapped, so this is purely a
  performance axis

### Against what it replaces

| build | median `t32` | behaviour |
|---|---:|---|
| pre-fix `ed59850` | 11.51 s | deterministic @ 4 |
| `f36def1` (hard-coded floor 2) | 9.31 s | deterministic @ 2 |
| **`e78757c` (bracket)** | **9.68 s** | **bimodal, 1 or 4** |

Reference points from the dense sweep: pin-1 8.61 s, pin-2 9.22 s, pin-3 8.98 s,
pin-4 10.42 s.

### The cost asymmetry makes it worse

| first decision | points | samples | bracket wall | ends at |
|---|---:|---:|---:|---|
| model accepted | 1 | 2 | 0.70–1.00 s | 1 ✅ |
| model rejected | 2 | 4 | **3.60–3.80 s** | 4 ❌ |

The failing path consumes the entire 4 s budget — roughly **36% of a ~10 s
run** — and then settles on the wrong split. The trajectory is
`4 → 1 → 4 → 5 → 4`: test the model answer, judge it worse, restore the opening,
test the adjacent candidate away from the model, fail to retain it, exhaust the
budget at the opening.

---

## 3. Mechanism

`OpeningBracketConfig::default()` uses `horizon: 200 ms`, which yields **2
samples** of evidence for the first decision.

At `t32` the true gap between producer 1 and producer 4 is 8.6 s vs 10.4 s —
about 20% — against run-to-run noise of roughly 7% (±0.6 s on a 9 s run). Two
samples cannot separate those reliably. An inconclusive first test is then
treated as a **rejection**, which triggers the expensive fallback search and
restores the opening.

`t8` and `t64` are stable because their gaps are large enough to clear the noise
at 200 ms; `t32` sits right at the boundary.

This is the same pattern flagged in the review of `f36def1`'s
`unresolved_floor_validation` — *inconclusive → revert → back to the opening is
the old bug gated behind noise*. It has now materialised as measured instability
rather than a hypothetical.

---

## 4. What is needed

1. **Require the same statistical separation to reject the model's answer that
   is already required to retain an alternative** (>5% plus interval
   separation). An inconclusive first test should **keep** the model's answer,
   not revert to the opening. The model is the primary signal; the opening is a
   guess with no evidential status. This is a small change and it targets the
   failure directly.

2. **Consider an asymmetric first horizon.** The first point gates a 0.8 s
   versus 3.8 s cost difference, so it is worth more evidence than subsequent
   points. Lengthening only the first horizon costs nothing in the common
   accept case if the accept decision is reached sooner.

3. **Re-gate `t32` on at least 8 repetitions**, requiring:
   - every run to land in producer 1–3 (no run at 4);
   - bracket wall time bounded on every run, not just the accepting ones;
   - median within 5% of the best measured pinned point.

4. **Add to the acceptance criteria generally: report `n` and the spread, not a
   single run.** `n = 1` is exactly what hid this. A bimodal outcome is
   invisible to a single measurement and this mechanism is, by construction,
   capable of bimodality — it makes a discrete decision from noisy evidence.

### Reproduction

```bash
cargo build --release --features rapidgzip-cpu-accounting   # 481b2be0…

for r in $(seq 1 8); do
  target/release/piscem-rs map-scatac -i /scratch1/rob/atac_test/chr12_idx \
    -1 /scratch3/rob/tmp/tb-atac-10x-R1.fq.gz \
    -b /scratch3/rob/tmp/tb-atac-10x-R2.fq.gz \
    -2 /scratch3/rob/tmp/tb-atac-10x-R3.fq.gz \
    -o /scratch3/rob/tmp/tbv-atac-t32-r$r/out -t 32 --quiet --dict auto
done
# then read thread_broker.opening_bracket.outcome and final_producer_limit
# from each map_info.json
```

---

## 5. On the constant-count pushback

The response declines to compress the surviving scATAC values to two, arguing
the four have independent semantics — opening hint, reader batch geometry,
progress publication resolution, sparse steady policy.

**That argument is accepted.** The review's "≤2" was a proxy for "stop encoding
performance beliefs in constraints", and that has been done properly: none of the
four is a safety floor, only two reach the broker configuration, and each now
carries provenance. Collapsing semantically distinct values to hit a count would
recreate the same category error in a different form. No further action wanted
here.

---

## 6. Housekeeping

* `dev` is **2 commits ahead of `origin/dev` and unpushed** (`0338bdc`,
  `e78757c`).
* `scripts/` is now showing as **untracked** where it was previously ignored
  per-file. Worth checking before anything is committed, so machine-local
  harnesses do not land in a source commit by accident.
* The only uncommitted working-tree change otherwise is a `.gitignore` line for
  these review documents.

---

## 7. Evidence

| what | where |
|---|---|
| `t32` × 8 reps (the bimodality) | `/scratch3/rob/tmp/tbv-atac-t32-r{1..8}/` |
| `t8` × 2, `t64` × 6 | `/scratch3/rob/tmp/tbv-atac-t{8,64}-r*/` |
| Flex `t8` regression check | `/scratch3/rob/tmp/tbv-flex-t8-r{1,2}/` |
| `t32` dense low-end sweep (1, 2, 3) | `/scratch3/rob/tmp/tb-lowend-p{1,2,3}-r{1,2}/` |
| pre-fix and `f36def1` comparisons | `/scratch3/rob/tmp/thread-broker-oracle-fill/`, `/scratch3/rob/tmp/tb-fix-t{8,32}-r*/` |

Binary `481b2be0163aa15c9feed1e1ac4f8c31d9f0c309447cf70809e3fa78f7380126`, built
from `e78757c` with `--features rapidgzip-cpu-accounting`. All runs on the
2,000,000-read / 528,410-mapped scATAC fixture unless stated.
