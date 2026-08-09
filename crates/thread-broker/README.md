# thread-broker

Divides a fixed thread budget between a **producer** (typically a parallel
decompressor) and the **consumer** that does the real work, by measuring the
running pipeline rather than by guessing up front.

> **Status: active development.** The control law is stable and validated, but
> this is pre-1.0 and the API may change. Several tuning constants are honest
> observations rather than derivations — `deadband_threads` most of all, whose
> own documentation says so. See [Limitations](#limitations).

## The problem

A user asks for `-t 32`. Those 32 execution slots have to be split between
decompressing input and processing it, and there is no good way to pick the
split ahead of time. It depends on the compression ratio, on how much work each
record costs, on how many input files there are, and on the machine. All of
those vary per run, and some vary *within* a run.

Guessing wrongly is expensive in both directions. Measured on real read-mapping
workloads: where a parallel decoder can add concurrency it is worth up to
**4.2×**; where it cannot, dedicating slots to it loses **8–28%**, because a
decode slot can only decode while an inline decoder is work-conserving.

The usual fix is a hill climber — nudge the split, keep what helps. That fails
here for a structural reason. The natural signal, "is the consumer starving?",
is **one-sided**: the set of splits where starvation is zero is a half-line, not
a point. Every split past the optimum reads as equally fine, so a climber walks
until it runs out of threads. An earlier version of this crate did exactly
that — it converged on 47 of 64 threads for decoding and ran 44% slower than the
right answer.

## The design

**Solve for the split; don't search for it.** The approach follows DS2 (Kalavri
et al., OSDI'18): estimate each stage's *true* rate — the rate it would sustain
if nothing throttled it — and compute the allocation directly, so a wrong split
is corrected in one move instead of a walk.

The key step is that the true rates are never needed in absolute terms. Over a
window in which `X` items cross the pipeline, the producer spends `X · s_p`
thread-nanoseconds and the consumer `X · s_c`. The ratio of *busy time* is
therefore the ratio of per-item cost, and `X` cancels:

```
d* = N · busy_p / (busy_p + busy_c)
```

No shared unit of work is needed, and neither side has to know what the other
counts. Each side reports only `Work { busy_nanos, items }`.

Around that sits a small state machine:

| phase | what it does |
|---|---|
| **Warm-up** | discards startup. Every workload measured reported ~99% consumer starvation in its first 75 ms — with no warm-up, *every* input looks producer-bound and the first move is always wrong |
| **Survey** | smooths several windows, solves, and moves only if the answer is outside the deadband |
| **Blackout** | discards windows spanning a move; a resize does work and moves no items, biasing the very ratio the decision came from |
| **Ratify** | gathers throughput evidence and reverts on a *regression*. Absence of improvement does not revert — the surface is flat near the optimum, so "no better" is the expected reading for a correct move |
| **Steady** | settles, then sleeps until the solved target drifts persistently. Watching the *model* rather than throughput catches a consumer becoming cheaper, which a throughput-based detector sleeps through entirely |

**Failures are advisory.** A refused resize, an acknowledgement timeout or a
sampler panic must never fail the surrounding run or corrupt its output. The
broker reports the failure; the caller decides.

## Using it

Implement `Consumer` and `Producer` for your two sides, then:

```rust
let broker = ThreadBroker::builder(consumer, producer)
    .budget(32)
    .initial_producer_slots(4)
    .build()?;
```

Each side reports `Work { busy_nanos, items }` as it runs, and is asked to
resize; the shrinking side always acknowledges before the other grows.

## Primary knobs

Defaults are the validated ones; every value below is a builder method.

| knob | default | what it trades |
|---|---|---|
| `budget` | — | total slots to divide. The one required setting |
| `initial_producer_slots` | — | where to open. A bad opening costs only the time to detect it |
| `deadband_threads` | `1` | how close to the model's answer is close enough. Moves are not free — each rebuilds per-thread state on both sides |
| `sample_interval` | `100 ms` | decision cadence |
| `warmup` | `400 ms` | startup discarded. Must span ≥ 2 sample intervals |
| `smoothing_windows` | `3` | noise rejection vs. responsiveness. `1` disables smoothing |
| `blackout_samples` | `4` | windows discarded after a move. Must be ≥ `smoothing_windows` |
| `ratify_samples` | `10` | evidence gathered before a move is kept |
| `regression_tolerance` | `0.05` | throughput loss that counts as a regression rather than noise |
| `resurvey_distance` | `1` | model drift, in threads, before a settled split re-opens |
| `resurvey_persistence` | `800 ms` | how long that drift must persist. Distance alone is budget-sensitive: with distance as the only gate, one run re-opened 18 times at `-t 64` while settling once at `-t 32` |
| `steady_state_policy` | `Responsive` | or `FreezeAfterConvergence` / `FreezeAfterFullCalibration` to stop adapting once settled |
| `opening_policy` | `Fixed` | or `Bracket(..)` to confirm a disagreement between the opening and the first stable answer |
| `min_consumer_threads` | `1` | floor for the consumer |
| `min_producer_slots` | `1` | **aggregate**, not per input. Raise only if your producer couples its limit to a one-time decision |
| `steady_probe_interval` | `sample_interval` | slower monitoring after convergence |

`EngagementPolicy` is separate and answers a prior question — whether to run a
parallel producer *at all*. It defaults to requiring 8 threads per producer
stream, since below that an inline, work-conserving producer wins.

## Limitations

- **`deadband_threads` is not derived.** The measurements originally justifying
  it were taken through a path that did not enforce the thread budget, which made
  a sharply peaked surface look like a wide plateau. Two threads is defensible
  against the corrected sweeps but a single constant is unlikely to suit every
  budget; it matters most when the budget is small.
- **`ratify_samples` assumes a stationary workload.** A deliberately changing or
  pausing workload is not a statistically identifiable ratification population.
  Validate the default's power against real windows before relying on it.
- **Producer measurement quality varies.** If the producer cannot report its own
  busy time, the broker falls back to sampling, which is noisier; it warns when
  it does.

## References

- Kalavri et al., *Three Steps is All You Need: Fast, Accurate, Automatic Scaling
  Decisions for Distributed Streaming Dataflows* (DS2), OSDI'18
- Suleman et al., *Feedback-Directed Pipeline Parallelism* (FDP), PACT'10
- Apache Flink FLIP-271 (autoscaler) and FLINK-31215
