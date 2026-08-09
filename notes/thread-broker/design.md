# `thread-broker`: dividing one thread budget between decoding and mapping

*Working document. Not tracked in git — see `.gitignore`.*

**Purpose of this document.** It describes the final design of
`crates/thread-broker`, names every type and function that carries load, and
records the mistakes made on the way to it. It is written to be read critically:
§11 lists the places I believe are weakest, and §12 lists what would falsify the
whole approach. Companion documents: `plan-adaptive-scheduling.md` (validation
strategy), `docs/threading-and-decompression.md` (the shipped, user-facing
behaviour).

**Status.** Post-audit implementation in progress. Safety, authoritative budget
planning, clean-window evidence, epoch-scoped rejection, robust ratification,
uncertainty-derived hysteresis, bounded telemetry, and synthetic nonlinear/noise
gates are implemented. The producer sampler is phase-aware and has direct
accuracy and same-binary overhead controls. Optional producer thread-lifetime
CPU accounting is merged upstream and has passed its direct hot-path and
application-level overhead gates. The full real modality oracle matrix and the
absence of released equivalents for the two git dependencies remain release
blockers. `audit.md` defines the gates; ignored
`completion-ledger.md` records live evidence.

---

## 1. The problem

A run gets one number from the user, `-t`. It must cover both the threads that
decompress the input and the threads that map the decompressed reads. Write the
split as `d` producer slots and `N - d` consumer threads.

The split matters more than the budget. On a 150 M read-pair FASTQ archive at a
fixed 64 threads, sweeping `d` was worth **5.75×**; doubling `N` bought 9%.

The right answer is not a constant. Per-thread mapping consumption ranged from
0.064 GB/s against a 96.3 M-k-mer index to 0.43 GB/s against a 1.5 M-k-mer one —
a 6.7× spread across workloads this tool ships for. Every attempt to *predict*
the split from a pre-run probe landed 7–60% off depending on the workload.

The throughput surface is a **narrow peak**, steep on both sides (§7):

| budget | optimum | −4 slots | +4 slots | +12 slots |
|---:|---:|---:|---:|---:|
| 32 | 12 | +29% | +13% | +115% |
| 64 | 22 | +17% | ~+3% | +9% |

Being wrong in either direction costs real money, and the tolerance is tighter
at smaller budgets. A controller has to land close, not merely on the right side
of a cliff.

### 1.1 What the environment provides

Both requirements below were unavailable when this work started and were obtained
upstream during it. That is why the design can assume them.

- **The consumer resizes mid-run.** `paraseq::parallel::ThreadPool` gained
  `set_threads`, `share(ways)` and `total_live()`
  ([fork](https://github.com/rob-p/paraseq), branch `feat/dynamic-thread-pool`,
  PR open). Workers can be added *and* retired while the pipeline runs.
- **The producer resizes mid-run, admission is decoupled from the runtime
  limit, and executing time is cumulative.** `rapidgzip_core::DecoderPool`
  (current reviewed revision `7b943ba`) gives several decoders one shared
  aggregate budget with a mutable `set_worker_limit` and an optional exact
  cumulative executing-region counter. The
  decoupling matters: on earlier versions a low limit at open time committed a
  decoder irreversibly to its sequential backend, which shipped broken three
  times before it was understood.

Both remain git dependencies until suitable crate releases exist; rapidgzip's
required changes are merged upstream, while the paraseq change remains on its
feature branch.

---

## 2. The control law that failed, and why it was structural

The first implementation was a four-quadrant hill climber: sample the consumer's
idle fraction, sample the producer's pressure, and if the consumer is idle *and*
the producer says our limit is holding it back, move threads to the producer;
if the consumer is busy and the producer is satisfied, move them back. Accumulate
evidence, decay it on disagreement, act at a threshold.

On Flex at `-t 64` it settled at **47 decode / 17 mapping** and ran **44% slower
than the best split**, with the trace reading `idle 6.5% → Starved → GrowProducer`
over and over.

This is not a tuning failure, and it is the single most important thing in this
document.

Let `starve(d)` be the fraction of consumer time blocked on an empty buffer. It
falls as `d` rises, reaches zero at the balance point — and then **stays** zero
for every larger `d`. The target set

```
{ d : starve(d) = 0 }
```

is a **half-line, not a point**. A controller driving `starve → 0` has no
restoring force anywhere on it. Any bias — and two independently sampled signals
always have some — walks it to the boundary while the measured error stays
genuinely near zero the whole way. The controller was not failing to converge. It
converged correctly, onto a set whose elements are not equally good.

A second, milder defect compounds it: the idle fraction is **scale-invariant**.
Halving both stages' speeds leaves it unchanged, so it carries no information
about *how far* to move. The step size had to be guessed, and `grow_step:
Proportional { cap: 8 }` was doing the real work.

Two independent literature reviews reached the same diagnosis. DS2 (Kalavri et
al., OSDI'18) §2 classifies exactly this controller — threshold policies over
utilisation, queue depth or observed rates — and gives the consequence as
"incorrect provisioning, oscillations, and long convergence times".

Kuo/Lim/Meerkov (1996) give the formal repair: the signed indicator
`e(d) = block₁(d) − starve₂(d)` has an **isolated zero** where `starve` alone has
a half-line. We do not use that indicator directly, but it is why a one-sided
signal can never work.

---

## 3. The control law that replaced it

### 3.1 True rate, not observed rate

DS2's central distinction. The **observed** rate of a stage is work over wall
time; the **true** rate is work over the time the stage spent *actually working*:

```
TRUE_PROCESSING_RATE = numRecordsInPerSecond / (busyTimeMsPerSecond / 1000)
```

The observed rate of a starved stage tells you how starved it is. The true rate
tells you how fast it *can* go, and is a property of the work rather than of the
current allocation. A controller built on wall time is feeding its own actuation
back into its measurement.

### 3.2 The closed form

With `s_p` and `s_c` the per-item costs in thread-time, a budget `N` split into
`d` and `N − d` finishes both sides together when `d / s_p = (N − d) / s_c`:

```
    d* = N · s_p / (s_p + s_c)
```

FDP (Suleman et al., PACT'10) expresses the same quantity as `τᵢ = coresᵢ /
avg_exec_timeᵢ` and picks the LIMITER as `argmin τ`; equalising the `τ` is this
equation.

### 3.3 The simplification that made it practical

The formula as written needs `s_p` and `s_c` in a **common unit of work**, which
forces the two stages to count the same items. That is a genuine problem here:
the decoder counts decompressed bytes and the mapper counts reads, so reconciling
them means summing full FASTQ record bytes on the mapping side and hoping no
input is uncompressed.

It is unnecessary. Over a window in which the pipeline moved `X` items, the
producer spent `X · s_p` thread-nanoseconds and the consumer `X · s_c`, so

```
    busy_p / (busy_p + busy_c)  =  s_p / (s_p + s_c)      exactly
```

The `X` cancels. **The ratio of the two stages' busy times over one window is the
cost share directly.** The implemented law is

```
    d* = N · busy_p / (busy_p + busy_c)
```

with no shared unit and no per-record byte accounting. It also handles a case the
byte form gets wrong for free: a run mixing gzip and plain inputs. The producer
simply does no work for the plain ones, its busy time is proportionally smaller,
and the split follows. Under the byte form the two stages would count different
populations of bytes and the answer would be biased by the compressed fraction.

This substitution is mine rather than the literature's, and it is the step a
reviewer should attack hardest. Its assumptions are stated in §11.1.

### 3.4 What production systems do with such a model

Not trust it. Flink's autoscaler (FLIP-271) computes a target parallelism from
true processing rates and then guards it:

- **`observed-scalability`** — record `(parallelism, true rate)` pairs; never
  assume a stage scales better than it has been seen to.
- **FLINK-31215** — a non-scalable bottleneck upstream makes the model keep
  scaling a stage that cannot go faster. Flink back-propagates it; we cap (§4.6).
- **`detectIneffectiveScaleUp`** — compare the achieved increase against the
  **model-predicted** increase, and revert if it falls far short. It ships
  **disabled by default**, which is the right calibration of how much weight one
  throughput comparison deserves.

FDP §3.3 is a hill climber, but in composition it is a *check* on a model-based
decision, with revert-on-regression and a tabu set — not the decision itself.

### 3.5 Rejected alternatives

| approach | why not |
|---|---|
| SPSA | degenerate in one dimension — the simultaneous perturbation *is* the gradient estimate |
| Kiefer–Wolfowitz | needs a smooth stochastic gradient; ours is a narrow peak |
| Golden-section search | irreversible interval elimination; under noise one bad sample discards the optimum permanently |
| Bayesian optimisation | a run yields a few dozen windows; sample cost exceeds the budget |
| Extremum seeking | requires persistent dither; every move here rebuilds multi-MB per-thread state |
| Probabilistic bisection (Horstein 1963; Waeber/Frazier/Henderson 2013) | genuinely applicable, `O(e^{-rn})` — but it answers a question we can solve in closed form. **Held in reserve** if the cost model proves unmeasurable |
| OSUB / unimodal bandits (Combes & Proutière, ICML'14) | same |
| COMPASS, SEDA, Dhalion | threshold or queue-length driven; the same one-sidedness as §2 |
| TBB | no auto-tuning of a producer/consumer split |

---

## 4. The implementation

`crates/thread-broker` — 750 lines of `lib.rs` (contract, config, builder,
instrumentation) and 690 of `controller.rs` (the loop). It knows nothing about
gzip or mapping, creates no workers, and owns exactly one thread: the sampler.

### 4.1 The measurement contract — `Work`

`lib.rs`, `pub struct Work { busy_nanos: u64, items: u64 }`.

`busy_nanos` must be time spent *doing the work*, blocking excluded. This is the
load-bearing field and the only real requirement the crate places on a caller.

`items` is a progress measure used **only** to compare throughput at one split
against another. It is never divided by the other stage's count, so the two
stages may count entirely different things. The type documents this at length,
because treating the two `items` fields as a shared currency is exactly the
mistake §3.3 avoids, and it is an inviting one.

Both fields are cumulative and must be monotonic. Before taking a delta, the
controller detects a decrease in either field and terminates adaptation with a
structured `WorkCountersRegressed { side, previous, observed }` failure and the
partial report. This is advisory at the piscem boundary, so mapping continues at
the last safe split. The remaining `saturating_sub` is defensive arithmetic,
not reset handling; silently treating a reset as zero work would stall or
misprice the controller until the new counter caught the old value.

### 4.2 The two traits

```rust
pub trait Consumer: Send + Sync {
    fn set_threads(&self, n: usize);
    fn live_threads(&self) -> usize;   // diagnostic only
    fn work(&self) -> Work;
}

pub trait Producer: Send + Sync {
    fn set_limit(&self, n: usize);
    fn limit(&self) -> usize;
    fn pressure(&self) -> ProducerPressure;   // vetoes growth; never sizes
    fn work(&self) -> Work;
}
```

`live_threads` is deliberately *not* the target: the gap between requested and
running is what separates "the broker decided badly" from "the broker decided
well and the pool has not caught up". It appears in the trace log and nowhere
else.

`ProducerPressure` has four variants — `Starved`, `Satisfied`, `SourceBound`,
`Inelastic`. The last two cap growth (§4.6); `Starved` has the bounded shrink-veto
exception in §4.5. None sizes the split. The enum's own docs carry the §2
argument, so that anyone tempted to size the split from it reads the reason not
to first. `Starved` used to carry a
`wanted: usize` demand estimate; nothing reads it under a solving controller, and
a demand estimate with no consumer to keep it honest is a liability, so it was
removed.

### 4.3 Instrumenting the consumer — `BusyMeter` / `BatchTimer`

`BusyMeter` is the shared counter; `BusyMeter::timer()` hands out a `BatchTimer`
that publishes **every `DEFAULT_FLUSH_EVERY` (256) items**, not once per batch.

That granularity is not an optimisation. A counter published once per batch cannot
be sampled faster than a batch completes: at 16 384 records and ~5 µs per record a
batch takes ~80 ms against a 100 ms window, so windows routinely contain no
completed batch, read a delta of zero, and report **maximum starvation for a
thread running flat out**. Measured before the fix: a deliberately expensive
kernel reported 100.0% idle — the exact opposite of the truth.

The unit must also be fast enough for the record cost. The real scATAC fixture
publishes completed-record progress every 64 records while a decision is open:
at 256, 25 ms windows were quantized into zero/256-record bursts and a genuinely
faster nonlinear point could not clear its confidence interval. Increasing the
number of such quantized windows is not the remedy; it extends calibration while
leaving the observation process bursty. Finer publication fixes the observable
signal itself. This cadence is phase-aware. Responsive steady state returns to
256; freeze drops the controller-side consumer at convergence, after which new
batches publish only on drop. The cadence is loaded once per paraseq callback
rather than per record. Only the completed-item counter uses the fine cadence.
Busy-time clock reads and updates remain at 256 records; the response test needs
fine throughput, not a new clock reading at every progress publication.

The first formal run exposed contention in a shared progress atomic. scATAC now
gives every processor a cache-padded progress shard. Its single writer loads the
cumulative value once per callback and uses relaxed stores at the 64-record
cadence; the controller sums shards only at its sampling cadence. Generic
consumers retain the shared additive 256-record meter. This preserves the
required resolution without multiplying read-modify-write contention.

`PISCEM_SCATAC_PROGRESS_FLUSH_EVERY` is a same-binary validation hook, not a
user-facing tuning control. The local runner requires exact first-position
balance after an earlier randomize/reverse implementation produced an invalid
20/10 split. The final 30-pair, 2-million-record 64-versus-256 run was exactly
15/15 balanced and retained canonical output/counts in all 60 runs. Paired
fine/generic medians were 0.99933 mapping wall, 0.99895 process wall, and
0.99965 aggregate CPU. Position-stratified one-sided 95% upper ratios were
1.00419, 1.00434, and 1.00441, respectively, so the formal <=1% sub-gate passes.

`BatchTimer` publishes on `Drop`, so an early return still records what was done,
and **must be constructed inside the per-batch callback**. That placement is what
keeps blocking out of the measurement: the callback is entered only once a filled
buffer is in hand, so the wait for that buffer falls *between* two timers rather
than inside one. A timer hoisted to thread scope would charge every `fill` to the
consumer's busy time, which is precisely the wall-time measurement the crate
exists to avoid. This is documented on the type; the four call sites are in
`src/mapping/processors.rs` (bulk PE, bulk SE, scRNA, scATAC).

### 4.4 The loop — `Phase`

`controller.rs`, `ThreadBroker::run`. A five-state machine sampled every
`sample_interval` (100 ms) after `warmup` (400 ms).

**Warm-up is load-bearing.** Every workload measured reported ~99% consumer
starvation in its first 75 ms, purely because its threads were still starting.
Without it, *every* input looks producer-bound and the first move is always wrong.
`makes_no_decision_during_warm_up` pins it.

| phase | does | leaves for |
|---|---|---|
| `Survey` | accumulate clean `Costs`, solve, and jump if the answer clears the uncertainty-derived deadband and has not failed in this epoch | `Draining` or `Steady` |
| `Draining { side, to, started, next }` | wait for the shrinking side's actual occupancy acknowledgement; fail on timeout | grow the other side, then `next` |
| `Blackout { left, next }` | discard `blackout_samples` **clean, stable-flow** windows, then `Costs::clear` | `next` |
| `Ratify { from, target, baseline, rates }` | compare like-for-like usable throughput blocks with pooled uncertainty | `Survey` if kept/inconclusive, `Blackout → Steady` if confidently reverted |
| `Steady { drifted_for }` | keep solving; re-open only on a drift ≥ `resurvey_distance` that persists for `resurvey_persistence` wall time | `Survey` |

`Draining` separates pool convergence from blackout. The shrinking side must
report actual occupancy at or below the request before the other side is grown,
so no pre-acknowledgement resize data can reach the model. `Blackout` then exists
because the move itself *does work* — a
new consumer worker builds its thread state and allocates a multi-MB record set —
which inflates busy time without producing items, biasing the very ratio the next
decision reads. Only usable windows with comparable buffer progress decrement
blackout. `build()` enforces `blackout_samples >= smoothing_windows`.

`Ratify` reverts **only on a confident regression** beyond
`regression_tolerance`, never on absence of improvement. Baseline and achieved
block uncertainty combine in quadrature for a one-sided 95% comparison;
overlap is recorded as inconclusive and creates no rejection evidence. Near the
optimum the surface is locally flat — 22 and 23
slots of 64 differ by 3%, inside this comparison's own noise — so "no better" is
an ordinary reading for a *correct* move, and treating it as failure would send
the broker back to wherever it happened to start, the one place with no evidence
for it. This is Flink's `detectIneffectiveScaleUp` lesson.

`OpeningPolicy::Bracket` is the opt-in correction for allocation-dependent stage
costs that the one-point service-cost model cannot identify. It is keyed on a
stable disagreement between the model answer and the opening, so a workload
where they agree pays zero points, samples, and wall time. The already-measured
opening is the baseline. The model answer is the first candidate; if it is
clearly retained, the experiment ends at the model answer. A loss that first
appears at the short opening horizon is extended to the ordinary ratification
sample count, so the common accept path stays short while an apparent rejection
gets more evidence. A regressed or inconclusive opening comparison triggers at
most two adjacent points, first away from the model and then toward it if the
first point fails.

Crucially, the opening is a pivot for candidate placement, not the retained
baseline. Each adjacent point is compared with the higher measured opening/model
rate and must prove both non-overlapping intervals and a greater-than-5% gain.
This high-water bar prevents a temporarily depressed model epoch from making a
second mediocre point look like an improvement. If no point clears it—or the
total deadline expires—the model remains the fallback. This
distinguishes t8, where producer five directly and decisively beats model point
one, from t32, where nonstationary startup can make the earlier opening block
look better without proving that any alternative beats the model. It also
implements the review's central rule literally at the final decision:
inconclusive evidence never restores the opening.

Persistent empirical cap evidence can end the ambiguity before local
exploration. If slack or source evidence independently establishes a useful cap
at or below the model target, buying adjacent producer slots contradicts both
the cost model and observed admission. The model is retained even when short
startup throughput blocks disagree. This is the measured t32 discriminator:
its cap is one, whereas t8 has headroom above model target one and therefore
still tests producer five. The rule depends on runtime evidence, not modality or
budget constants.

The default bracket is three points, 200 ms of evidence per point, and a four-
second wall budget. These are generic bounds rather than scATAC constants. The
controller will not begin another point without the remaining evidence and
blackout horizon, wakes at the deadline, restores the last proven split on
expiry, and reports `points_measured`, `samples`, `wall_nanos`, and `outcome`.
An in-flight pool shrink must still obey its independent two-second resize
acknowledgement contract before the other side can grow. The bracket is never
rearmed by a steady-state resurvey, so it changes no settled cadence or recurring
work. `FreezeAfterConvergence` deliberately skips it;
`FreezeAfterFullCalibration` completes it before teardown.

The paired synthetic regression is intentionally not satisfiable by trusting a
single opening: both shapes open at producer four and produce model answer one,
but the mapping-heavy t32-like shape retains one while the negative-scaling
t8-like shape rejects one, samples producer five, and retains it.
This replaces the old geometric 4→6→7→5 search and its budget-specific scATAC
constants with the question the controller actually needs to answer.

After calibration, scATAC defaults responsive monitoring to a 5 s probe
interval because measured stable runs do not change their selected split. This
does not affect warm-up, response probing, ratification, or resurvey cadence;
`PISCEM_THREAD_BROKER_PROBE_INTERVAL_MS=25` is the explicit same-binary override
for applications that need normal-resolution regime-change response.

After a *kept* or inconclusive move the loop returns to `Survey`. DS2 reports
convergence in at most a few such rounds; uncertainty-derived hysteresis and
epoch-scoped rejection stop it cycling.

### 4.5 The solver — `Costs::solve`

`Costs` is a ring of `(consumer_busy, producer_busy)` pairs, `smoothing_windows`
deep. `solve` sums both — sums, not an average of per-window ratios, so a window
in which the pipeline barely moved contributes proportionally rather than equally
— and returns

```rust
let share  = producer as f64 / total as f64;
let ideal  = (budget as f64 * share).round().max(1.0) as usize;
let target = ideal.clamp(lo, hi).min(caps.useful().max(lo));
```

with `lo = min_producer_slots`, `hi = budget − min_consumer_threads`. It returns
`Solved { target, snapshot: Model { producer_cost_share,
producer_cost_share_uncertainty, ideal_producer_slots, useful_cap,
useful_cap_reason, effective_deadband_threads,
effective_resurvey_distance } }`.

`Model` is public and reaches the user through the CLI log line:

```
thread broker: decode is 37% of per-read cost -> wanted 12 slots, usable ceiling 13
```

That is deliberate. The cost share is a **directly falsifiable claim about the
workload**, and it separates "the model was misinformed" from "the model was
overruled by the cap" without a debugger. It is the first thing to look at when a
split looks wrong.

One bounded correction protects the constant-cost assumption without restoring
hill climbing. If the model asks to *remove* a producer slot while the producer
still reports runnable work queued behind its current limit (`Starved`), that
shrink is vetoed and counted in the report. Pressure never chooses a larger
target. This matters when one slot under-reports the thread-time or parallelism
that a second slot unlocks; the real `-t 8` short-run surface measured a 22%
throughput gain at two slots despite a one-slot model share near 5%.

### 4.6 The empirical cap — `Caps`

The model assumes a stage can use every thread it is given. A producer reading
from a saturated disk cannot, and the cost share cannot tell — it reports how much
*work* decoding is, not whether there is anything to decode from. Uncapped, the
solver buys slots that sit idle, which is FLINK-31215.

`Caps::observe_slack` records `(achieved_concurrency, granted_limit)` only for
stable-flow windows, where achieved concurrency is `dp.busy_nanos / window`.
`Caps::observe_pressure` independently records direct source/elasticity pressure
whenever the producer is not resizing; downstream buffer drift does not make
that admission classification less true. `Caps::useful()` returns the
minimum of two one-sided bounds — a cap is only ever recorded when the producer
*demonstrably* failed to use what it had:

1. **Slack.** More than a thread's worth of unused grant must persist for
   `cap_persistence` wall time before it establishes the best achieved
   concurrency plus one.
2. **Pressure.** `SourceBound` or `Inelastic` must persist for the same duration
   before capping at the current limit.

The cap is duration-based, not sample-count-based. Evidence expires after
`cap_history` (800 ms by default) of contradictory observations, so a 0.25-1 s
lull cannot hold down later decode demand for more than one second. `Model`
reports `slack`, `source`, or `slack_and_source` as the reason.

### 4.7 Epoch-scoped rejection evidence

A confidently reverted target is recorded with its workload epoch, baseline and
achieved rates and uncertainty, and the cost share that proposed it. The same
target is not retried while that model remains unchanged; the last safe split is
kept instead of walking through unsupported neighbours.

A persistent cost-share or throughput regime change increments the epoch and
clears old rejection evidence. This preserves termination within one regime
without making the broker permanently deaf when a later regime legitimately has
the same optimum. The deterministic two-regime test revisits a rejected target
within the five-moves-per-epoch gate.

### 4.8 Change detection — drift, not CUSUM

The plan called for a CUSUM on throughput: settle, then sleep until surprised. It
was implemented, and a test killed it.

A throughput CUSUM is **one-sided**. When the consumer suddenly gets *cheaper* — a
new file, a different sample, a smaller index — throughput does not fall at all:
the pipeline stays pinned at the old producer-limited rate while a much better
split goes unclaimed. The monitor sleeps through it.
`re_solves_after_a_regime_change` fails against that implementation.

Watching the **solved target** drift detects both directions with one mechanism,
because the target is computed from the cost ratio rather than from the outcome.
Hysteresis is derived from the observed 95% cost-share uncertainty at the current
budget, with one-thread absolute floors. Re-opening additionally requires 800 ms
of consecutive drift evidence, making it harder than an initial move without making
a one-slot regime change invisible at small budgets. `build()` rejects a
resurvey floor below the movement floor, since a settled split that immediately
re-opens itself oscillates forever.

**Distance alone was not enough, and only a real run showed it.** Noise in the
solved target scales with the budget; a fixed distance does not:

| budget | moves | resurveys | converged |
|---:|---:|---:|---|
| 32 | 1 | 0 | yes, 1.4 s |
| 64 | 16 | 18 | **never** |

Identical workload, identical noise in the *share*, twice the noise in threads.
`resurvey_persistence` (800 ms by default) requires the drift to **persist** in
wall time, so changing the observation cadence does not silently redefine the
guard. The first sparse observation earns only one active-cadence window; it
cannot claim the entire sleep since the preceding clean probe as corroboration.
After the fix, `-t 64` ran 3 moves, 0 resurveys, converged in 2.6 s.

When a nonlinear point is retained, an upside-only common-mode throughput
change at that same split does not reopen the curve: the scATAC consumer set
warms after shrinking and otherwise reopened repeatedly without evidence that
another split was better. A confident loss still reopens. Epoch-scoped model
rejections retain two-sided rate change as an evidence-expiry mechanism.

### 4.9 Window rejection — `usable_window`

A cost share is a ratio, and a window containing a handful of nanoseconds of work
gives a ratio dominated by where the boundaries fell relative to the counters'
flush granularity. The CLR thread pool's hill climber guards the same hazard from
the other side, rejecting a window when in-flight work is large relative to
completed work.

```rust
fn usable_window(dc: Work, dp: Work, window: Duration) -> bool {
    let floor = window.as_nanos() as u64 / 20;
    dc.items > 0 && dc.busy_nanos.saturating_add(dp.busy_nanos) > floor
}
```

The test is on the **total**, not on each side. That is deliberate: a producer
doing no work at all must be admitted, not discarded. A run whose inputs are
already uncompressed has no decode cost, and the honest reading is a share of
zero, which solves to the producer floor and hands the whole budget to the
consumer (`hands_the_whole_budget_over_when_there_is_nothing_to_decode`).
Requiring both sides to be busy would reject every window forever and strand the
initial guess.

Progress is required as well as busy time for the **cost ratio**: a window in
which the pipeline moved nothing says only that something is stalled, and which
side happened to be spinning during it is not evidence about relative costs.
The same zero is nevertheless essential **throughput** evidence. Dropping it
conditions a bursty rate on observing a completion and systematically inflates
slow splits, so rate histories and ratification retain zero-progress windows.

### 4.10 Applying a split — `begin_move` / `Draining`

Always requests shrink, observes actual occupancy acknowledgement, and only then
grows. Growing first, or merely ordering two setters, lets asynchronous retirement
briefly sum above the budget. Refusal and timeout are surfaced as typed errors
with a partial terminal report.

The direction comes from the controller's owned `(from, to)` split rather than from
`Producer::limit()`. Reading it back would reintroduce the confusion the broker
avoids by tracking its own view of the split: a pool that has not converged to the
last request reports the old value, and the ordering decision would be made
against a split that no longer exists.

### 4.11 Configuration

`BrokerConfig` with a hand-rolled builder (`ThreadBroker::builder` /
`builder_with`). Defaults: 100 ms interval, 400 ms warm-up, 3 smoothing windows,
deadband floor 1, blackout 4, ratify 10, fixed opening policy (an opted-in
bracket defaults to three points, 200 ms per point and four seconds total),
regression tolerance 0.05, resurvey
distance floor 1, resurvey samples 8, 2 s resize timeout, 800 ms cap history,
300 ms cap persistence, bounded trace capacity 256, and floors of 1 and 1.
Piscem uses 25 ms / 100 ms for effective budgets up to eight so three clean
windows fit before a quarter-second job ends; larger budgets use the defaults.
The small-budget adaptive opening is two slots (`N/4`) rather than one (`N/8`),
with the pressure veto above preventing the model from discarding demonstrated
runnable parallel work. Mapping-heavy input remains free to reclaim the slot.
scRNA and Flex use a measured quarter-budget opening. scATAC uses one four-slot
opening at every budget and opts into `OpeningPolicy::Bracket`; its producer
safety floor remains the shared-pool value one. This removes the `budget <= 8`
threshold and lets the same mechanism select producer five at t8 and producer
one at t32/t64. Fixed and serial requests are unchanged.

Steady-state behavior is a separate builder choice:

- `SteadyStatePolicy::Responsive` (the default) keeps the regime-change logic.
  `steady_probe_interval(Duration)` changes only the cadence after one ordinary-
  cadence clean steady window has established convergence. Long waits are
  interruptible, so `finish()` does not wait for the next probe.
- `SteadyStatePolicy::FreezeAfterConvergence` takes the model-only path: it pays
  warm-up, survey, guarded model moves, blackout, ratification, and clean
  convergence, but deliberately skips opening calibration. It
  then terminates the controller and releases both recurring producer sampling
  and fine consumer progress publication, at the explicit cost of missing
  negative scaling and later regime changes.
- `SteadyStatePolicy::FreezeAfterFullCalibration` first performs responsive
  mode's bounded opening bracket, then stops the controller,
  producer sampler, and fine publication at clean convergence. It retains the
  zero-recurring-cost property while paying the full startup search needed by
  response curves such as scATAC's.

Piscem exposes
`PISCEM_THREAD_BROKER_POLICY=responsive|freeze-after-convergence|freeze-after-full-calibration`
and `PISCEM_THREAD_BROKER_PROBE_INTERVAL_MS=N` as same-binary validation hooks.
Applications should normally choose through the builder. Reports include the
policy, effective probe interval, controller lifetime/sample count, and whether
monitoring stopped at convergence.

`build()` rejects six inconsistent combinations, each of which is a bug someone
would otherwise ship:

- budget < 2, zero sample interval;
- `warmup < 2 × sample_interval` — a broker that decides on one sample decides on
  the startup transient;
- `blackout_samples < smoothing_windows` — the next decision still averages
  pre-move windows;
- `resurvey_distance < deadband_threads`;
- invalid cap durations, buffer drift fraction, resize timeout, or trace size;
- floors that do not fit in the budget.

The earlier API had two footguns, both fixed:

- A `.config()` method that **replaced** everything set before it, silently
  wiping a floor of 4 back to 1. Replaced by `builder_with(consumer, producer,
  config)`, which takes the config first so there is nothing to wipe.
- Defaults that **failed the crate's own validation** (250 ms interval against a
  150 ms warm-up). Every unit test passed an explicit config, so the defaults were
  the one combination never exercised, and the call site swallowed the error with
  `.ok()`. Now there is `the_default_configuration_builds`, and construction and
  broker thread completion return errors to the mapping command.

---

## 5. The adapters

The crate is generic; everything workload-specific is in `src/io/broker.rs`.

**`MappingConsumer`** wraps `paraseq::parallel::ThreadPool` and holds an
`Arc<BusyMeter>` cloned from `MappingStats` — owned rather than borrowed, because
the broker runs on its own thread for the length of the job.

`live_threads` returns `pool.total_live()`, the **aggregate**, not this handle's
share. `Collection` splits the pool across readers and the parent's own share never
runs anything, so `live()` reads zero for the whole run and every consumer measures
as 0% idle however starved it is. That was a real bug; the fix is upstream in the
fork at `40ae74f`.

**`DecodeProducer`** wraps `rapidgzip_core::DecoderPool` plus its handles.

Each decoder handle publishes an exact cumulative executing-region integral.
The production adapter sums those lock-free counters at the controller's
existing observation cadence; it creates no polling thread and performs no
high-frequency clock read in piscem. The rapidgzip feature adds a monotonic
clock read, wrapping live-time balance updates, and a packed transition epoch
at each already-existing decoder-task CPU-region begin/end boundary. The epoch's
low half counts boundaries in progress and its high half counts completed
boundaries. A reader accepts a live snapshot only when the epoch is unchanged
with no transition in progress and both balance/busy values are independently
stable and plausible. This prevents both torn updates and a complete
cross-counter ABA, either of which could otherwise make a cumulative result
temporarily overlarge and later regress. Its feature-off build compiles all of
these operations and branches out. Thirty paired direct-decoder comparisons
passed the no-measurable-overhead gate with byte-identical output: CPU ratio
median/one-sided-95%-upper `0.995767/0.998581`, wall
`0.998233/1.000000`.

A completed-region-only counter was also implemented and rejected. Although
monotonic and cheap, its publication lag caused 17--313% p95 error at 25--100 ms
on real sparse, bursty, and stored paths. Exact live accounting is necessary for
the decision window; thread-lifetime counters remain suitable only for
whole-run reconciliation.

The earlier phase-aware `BusySampler` remains only as a compatibility and
validation path. It sums instantaneous lock-free `busy_workers` counts and
integrates endpoints with a trapezoid: jittered around 3 ms during calibration,
25 ms in monitoring, and roughly four observations per explicitly sparse
steady probe. It deliberately does not use pool-permit occupancy, because a
reusable task may retain a permit across a nonblocking handoff while its
executing flag is clear. Sampled-versus-exact real-archetype validation showed
why it is not the production signal: sparse and bursty decision windows can
alias even when whole-run error is small.

For producer component accounting, rapidgzip now has a compile-time optional
`cpu-accounting` feature. It reads the calling thread's CPU clock at thread
registration and thread exit and performs one relaxed cumulative update at
exit. It adds no clock read, counter update, or branch to decode task begin/end
or inflate loops. These lifetime counters deliberately are not substituted for
the executing-worker signal used by the controller: they include all CPU work
on worker and auxiliary threads and omit still-running threads, making them a
whole-run component/reconciliation signal rather than a decision-window signal.

`pressure()` is the judgement-heavy part, because **upstream deliberately does not
report a source-bound state** — its own docs say it cannot yet separate source
waiting from inflate time reliably enough to promise as a contract. So
`SourceBound` is *inferred*, and the direction of the inference is chosen from the
asymmetry of the two errors:

- calling a genuinely starved decoder source-bound makes the broker reclaim decode
  threads — worth 22–48% when decode is really the constraint;
- calling a genuinely source-bound decoder starved wastes a few threads that will
  idle — a few percent of CPU.

So source-bound is claimed only on strong evidence (nothing queued, nothing
executing, our limit demonstrably not the constraint), and everything short of
that resolves toward `Starved`.

`tests/broker_adapters.rs` exercises this against real decoder states —
`producer_work_is_bounded_busy_time` checks the integrator accumulates at all and
never exceeds `limit × elapsed`, since a counter stuck at zero would silently
price decoding as free and hand the entire budget to mapping on every workload.

The same two adapters are what salmon, cuttlefish 3 and ruSTAR would write. The
control law should need nothing.

---

## 6. The tests

`crates/thread-broker/tests/control_law.rs`, 16 tests against a fake whose optimum
is known in closed form. Every assertion is against `Pipeline::optimum(budget)`
computed from the fake's own parameters, never a remembered number.

**Two earlier versions of this file tested nothing**, which is worth recording
because both looked reasonable:

- the first had a background thread sleeping 200 µs per item, which Linux will not
  honour;
- the second let the consumer's idleness be set independently of the producer's
  capacity, so the broker could move every thread it liked and the "starvation" it
  was chasing never responded.

Both passed against a controller that was badly wrong. What makes the current one
a test is that supply and demand are **joined**:

```
supply     = producer_slots / s_p     (optionally capped by the source)
demand     = consumer_threads / s_c
throughput = min(supply, demand)
```

and *both* busy counters are derived from that one throughput, so the broker
cannot improve one reading without paying for it in the other. The simulation is
advanced by polling (`Pipeline::advance`) rather than by a thread, so it stays
exact however the scheduler treats the test process.

| test | what it pins |
|---|---|
| `solves_a_{balanced,decode_heavy,mapping_heavy}_split` | the closed form at three cost ratios |
| `does_not_walk_away_from_the_optimum` | **the §2 regression.** Start *at* 32/64 and stay; the old law drifted to 47 from exactly here |
| `reaches_a_distant_optimum_in_a_handful_of_moves` | 2 → 48 of 64 in ≤3 moves; distance must not set the cost |
| `does_not_buy_slots_the_source_cannot_feed` | §4.6, FLINK-31215 in miniature |
| `reverts_a_move_that_makes_throughput_worse` | the model made wrong on purpose (20× producer bias); only throughput can catch it |
| `settles_and_stops_moving` | the tail of the trajectory is flat |
| `re_solves_after_a_regime_change` | §4.8; **fails against a throughput CUSUM** |
| `hands_the_whole_budget_over_when_there_is_nothing_to_decode` | §4.9 |
| `never_breaches_either_floor`, `never_oversubscribes_the_budget` | the contract with the caller |
| `makes_no_decision_during_warm_up` | §4.4 |
| `the_default_configuration_builds`, `rejects_a_blackout_shorter_than_the_smoothing_window` | §4.11 |
| `the_meter_publishes_within_a_batch` | §4.3 |

---

## 7. Measured result

Full 10x Flex v2 dataset, human index, `map-scrna`, AMD EPYC 9575F. Decode slots
are the *total* across both mates. Mapped-read counts were identical
(1 994 961 121) at all 13 points, so this is purely a performance axis.

**`-t 32`:**

| decode slots | wall | avg threads used |
|---:|---:|---:|
| 8 | 235.27 s | 22.0 |
| **11 — `auto` chose this** | **191.70 s** | 28.2 |
| **12 — best in sweep** | **182.52 s** | 29.6 |
| 16 | 206.41 s | 25.7 |
| 20 | 266.54 s | 19.6 |
| 24 | 393.32 s | 13.1 |

**`-t 64`:**

| decode slots | wall | avg threads used |
|---:|---:|---:|
| 16 | 132.86 s | 47.5 |
| **22 — best in sweep** | **113.08 s** | 55.3 |
| **23 — `auto` chose this** | **116.89 s** | 55.0 |
| 28 | 118.49 s | 51.9 |
| 32 | 123.28 s | 47.6 |
| 40 | 151.31 s | 36.1 |

| budget | `auto` | best | gap | moves | reverts | resurveys | settled after |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 32 | 11 | 12 | **5.0%** | 1 | 0 | 0 | 1.4 s |
| 64 | 23 | 22 | **3.4%** | 3 | 1 | 0 | 2.6 s |

The control law it replaced, same workload: **44% off, 16 moves, never converged.**

Two observations. The optimum coincides with peak machine utilisation (29.6 of 32),
and utilisation falls away monotonically in both directions — the signature of a
budget genuinely divided rather than spent twice. And the `-t 64` row is the direct
test of §4.8: the same build before the persistence fix ran 16 moves and 18
resurveys and never converged.

The measured cost share was 37% at `-t 32` and 26% at `-t 64`. The original
interpretation that these must agree was too strong; §11.2 now separates
same-split instrumentation repeatability from real allocation-dependent cost.

---

## 8. Pitfalls, in the order they were hit

Each of these cost real time; several passed a test suite first.

| # | pitfall | how the design answers it |
|---|---|---|
| 1 | **Wall time as the signal.** Idle fraction is endogenous in the split being controlled, and its zero set is a half-line. | Busy time only, and a solved split rather than a searched one (§2, §3). Enforced by the `Work` contract. |
| 2 | **Per-batch publication.** A counter updated once per batch reads zero in windows containing no batch boundary — maximum starvation for a thread running flat out. Measured: 100.0% idle for the busiest kernel. | `BatchTimer` publishes every 256 items (§4.3). `the_meter_publishes_within_a_batch`. |
| 3 | **The startup transient.** Every workload read ~99% starved in its first 75 ms; without a warm-up, every first move is wrong. | Explicit `warmup`, validated to span ≥2 intervals. `makes_no_decision_during_warm_up`. |
| 4 | **Reading the split back from the pools.** Conflates "what I asked for" with "what has taken effect", so two consecutive decisions double-spend the same threads. | The broker owns `split` and never reads it back — including in `apply`, where only the *direction* was still read back until this cleanup (§4.10). |
| 5 | **`total_live()` vs `live()`.** `Collection` shares split the pool; the parent's share runs nothing, so the consumer measured 0% idle for the whole run however starved. | `MappingConsumer::live_threads` uses the aggregate. Fixed upstream at `40ae74f`. |
| 6 | **Evidence that never decayed.** A balanced run still throws occasional starved samples; with nothing opposing them the split drifted 30 → 4 consumer threads. | Obsolete under a solver, but the general lesson — a controller must have a restoring force in *both* directions — is why §4.8 watches drift rather than a one-sided CUSUM. |
| 7 | **Deciding per sample.** Jitter clusters, so two consecutive dips reach any threshold. | `smoothing_windows`, and sums rather than averaged ratios (§4.5). |
| 8 | **A `.config()` builder method that wiped earlier settings**, silently resetting a floor of 4 to 1. | `builder_with(consumer, producer, config)` takes the config first (§4.11). |
| 9 | **Defaults that failed the crate's own validation**, swallowed by `.ok()` at the call site. Every test passed an explicit config. | `the_default_configuration_builds`; call sites use `.inspect_err`. |
| 10 | **Test fakes that could not fail.** A 200 µs sleep Linux will not honour; then an idleness signal independent of producer capacity. | The coupled `Pipeline` (§6). |
| 11 | **A throughput CUSUM is one-sided** and sleeps through a change that creates upside. | Watch the solved target drift (§4.8). |
| 12 | **A fixed drift distance does not scale with the budget** — 18 resurveys at `-t 64`, 0 at `-t 32`. | Uncertainty-scaled distance plus wall-time `resurvey_persistence` (§4.8). |
| 13 | **A running-minimum cap never recovers**; one lull pins decode concurrency for the run. | Persistent, duration-based cap evidence that expires after contradictory observations (§4.6). |
| 14 | **An asymmetric ratify comparison.** The baseline was one pre-move window against `ratify_samples` post-move windows — one unlucky window sets a bar the move cannot clear. | `Recent`, averaging the same span on both sides (§4.4). Found in this cleanup pass. |
| 15 | **`--decoder parallel=N` never reduced `map_threads`**, so a pinned request doubled the budget instead of dividing it: 41.3 average threads against a budget of 32. | Fixed in all three CLIs (§9). |
| 16 | **Believing an "oracle" that used more machine.** The 41.3-thread run looked 41% faster than `auto` and nearly buried a correct model. | §9's rule: divide CPU seconds by wall time, check against the budget, *then* compare. |
| 17 | **`claim_first_slot` was unreachable**; mutation testing showed no test could distinguish it. | Deleted rather than kept and tested. |

---

## 9. Two invalidated measurements

Both are recorded because the *conclusions drawn from them* shaped the design.

**The pinned-path oversubscription bug.** `--decoder parallel=N` overwrote
`plan.decode_budget` without touching `plan.map_threads`, so both sides were spent
in full. At `-t 32`, `--decoder parallel=16` ran an average of **41.3 concurrent
threads** against a budget of 32 and looked 41% faster than `auto` purely by using
half a machine again. This is a shipped-behaviour bug — a user who pins the decoder
gets a run that ignores `-t` — fixed in all three CLIs by recomputing
`map_threads = num_threads − decode_budget`. It deserves a release note.

**The "flat basin" was mostly an artefact.** An earlier `-t 64` sweep recorded
everything from 28 to 60 decode slots landing within 8%, and several decisions
leaned on it, including the framing that "precision near the optimum is nearly
worthless". That sweep went through the same unenforced path, so its high-decode
points never paid for their decode slots and a peaked surface was flattened into a
plateau. Corrected: 40 slots is **34%** worse than 22, not 8%.

Part of the effect is real — the peak *is* broader at 64 than at 32 (within 8% for
roughly 22–30 slots, against 11–13 at `-t 32`), because a larger budget has more
absolute slack either side of the balance point. What was not real is the extent.

Consequences: the old table is void; `deadband_threads` rests on reasoning now
known to be wrong even though the value survives (§11.3); and the §7 result is a
stronger claim than it looks, because the neighbouring sample points are 13–29%
worse rather than within 8%.

**The rule, now three for three** (this, the `PISCEM_RAPIDGZIP_THREADS` mistake
where the env var set only the *starting* split while the broker overrode it, and
the index-load confound that inverted two conclusions): **divide CPU seconds by
wall time and check the result against the budget before believing any
comparison.**

---

## 10. Where the code lives

| what | where |
|---|---|
| `Work`, `Consumer`, `Producer`, `ProducerPressure` | `crates/thread-broker/src/lib.rs` |
| `BrokerConfig`, `ThreadBrokerBuilder`, `build()` validation | same |
| `BusyMeter`, `BatchTimer`, `DEFAULT_FLUSH_EVERY` | same |
| `ThreadBroker::start` / `run` / `begin_move`, `Phase` | `crates/thread-broker/src/controller.rs` |
| `Costs::solve`, `Solved`, `Model` | same |
| `Caps::observe_slack` / `observe_pressure` / `useful`, duration-based cap evidence | same |
| `Recent`, `usable_window`, `rate_per_second` | same |
| `BrokerReport`, `RunningBroker` | same |
| control-law tests, the coupled `Pipeline` fake | `crates/thread-broker/tests/control_law.rs` |
| `MappingConsumer`, `DecodeProducer`, `BusySampler` | `src/io/broker.rs` |
| adapter tests against a real decoder | `tests/broker_adapters.rs` |
| `BatchTimer` call sites (4 batch loops) | `src/mapping/processors.rs` |
| broker construction, report logging, pinned-path budget | `src/cli/map_{bulk,scrna,scatac}.rs` |
| `plan_thread_budget`, `open_input_pooled` | `src/io/fastx.rs` |
| decoder preference, FIFO/non-rewindable guard | `src/io/calibrate.rs` |

---

## 11. Remaining weaknesses after the audit

Listed for a reviewer to attack. I believe each is a real exposure, not a
hypothetical.

### 11.1 The busy-ratio substitution assumes steady state

`busy_p / (busy_p + busy_c) = s_p / (s_p + s_c)` holds exactly when both stages
processed the **same `X` items** in the window. Over a window where the buffer is
filling or draining they did not, and the estimate is biased by the difference.

The producer now exposes aggregate buffered progress. Material fill/drain windows
are rejected, and blackout/ratification advance only on stable-flow evidence.
The deterministic transient test proves the broker makes no move from mismatched
work. The final remote-pinned real alternating/source-bound Gate G matrix passed
all six 25/100 ms cells with only 0.0003--0.0029 percentage points of exact
whole-run share bias, persistent-cap identification in 0.376--0.400 s, canonical
output equality, and the eight-slot budget invariant.

### 11.2 Producer measurement and allocation-dependent service cost

Measured before the audit: **37% at `-t 32`, 26% at `-t 64`** on the same file.
The original document incorrectly required budget invariance. Under cache,
memory-bandwidth, and serialized-work limits, per-record thread-time may really
depend on allocation. Repeatability at the same pinned split and agreement with
cumulative ground truth are the correct instrumentation checks.

Candidate explanations, none yet tested:

- **The producer's busy time undercounts.** `busy_workers` counts pool slots;
  the decoder's coordinator and scanner threads (`auxiliary_threads`) do real CPU
  work in no slot at all, as does the copy into the reader's buffer. If that
  uncounted fraction is roughly constant per unit of *data*, it inflates the share
  more at small budgets, where the counted part is smaller. Direction plausible,
  magnitude unknown.
- **The consumer's busy time overcounts under contention.** `BatchTimer` measures
  wall time between publishes, so a descheduled thread is charged as busy. At
  `-t 64` on a 128-core box that should be negligible; at `-t 32` more so. This
  predicts the *opposite* sign, so it is probably not the explanation.
- **Sampled-counter bias is closed for production.** The controller now reads
  rapidgzip's exact cumulative executing-region counter. The compatibility
  sampler remains measurable but is not constructed on the production path.

**This remains the most important open question.** Gate F pins at least five
splits per workload and records both stage costs, throughput, CPU, and buffer
occupancy. If cost varies by more than 10%, held-out throughput predictions must
still be within 10%. The synthetic fake now makes both stages deliberately
sublinear and `auto` reaches at least 90% of its discrete oracle in five moves;
real modalities remain unmeasured.

### 11.3 Hysteresis now derives from uncertainty

The model retains a 95% half-width across clean window shares. Effective movement
and reopen thresholds scale that uncertainty by the budget, with one-thread
floors. Across budgets `{2,4,8,16,32,64,96}`, 1,000 seeded 30-second stable
traces have the required >=99% no-false-resurvey rate, and >=99% of feasible 10%
budget shifts reopen within two seconds. The remaining caveat is calibrating the
synthetic noise level from real window variance.

### 11.4 Rejection is epoch-scoped

The permanent absolute set is gone. Rejections carry evidence and expire only
after a persistent model/throughput change opens a new epoch. The two-regime
revisit test passes; more than two epoch transitions and isolated zero-progress
ratification windows remain useful adversarial additions.

### 11.5 `SourceBound` is inferred, not reported

Upstream declines to promise the distinction (§5). The inference remains
heuristic, but cap evidence now requires 300 ms persistence, records its reason,
and expires after 800 ms of contradiction. Deterministic 25/100 ms cadence tests
identify a true ceiling within one second and recover from a one-second lull
within another second. The real source-bound fixture is still required.

### 11.6 One workload, two thread counts

bulk SE, bulk PE, scRNA (non-Flex) and scATAC are unvalidated, as is `-t 8`. Flex
is the *most* decode-bound workload we ship; a mapping-bound one exercises the
opposite side of the surface, where the correct answer is a producer floor of 1–2
and the risk is over-buying decode. The design is not validated there.

### 11.7 The producer's `items` is bytes and the consumer's is records

Correct by the contract (§4.1) — they are never divided by each other — but the
crate cannot *enforce* it, and an implementor who assumes a shared currency gets a
silently wrong split scaled by the mean record size. The only defence is the
documentation on `Work`. A newtype per stage would make it unrepresentable.

### 11.8 Sampler and optional upstream accounting overhead are measured independently

`PISCEM_FIXED_DECODE_MEASUREMENT=off|monitoring|calibration` is an
environment-only validation hook for an aggregate-pinned run. All three modes
use the same binary and exact fixed mapping/decode split; only the sampler
changes. The local overhead runner uses counterbalanced randomized blocks,
canonical output digests, process CPU, and a one-sided 95% bootstrap equivalence
test. Separately, rapidgzip's compile-time optional lifetime CPU accounting is
tested with two independently compiled binaries over paired, counterbalanced
direct-decoder runs. Its feature-off build compiles the clock reads and counters
out entirely. Both the direct decoder and piscem feature-on/off comparisons use
a one-sided 95% upper bound of <=1% for wall and CPU time. Raw results and the
completion ledger record the measured bounds.

### 11.9 Whole-broker cost must be fixed and fractional

Sampler-only cost is not the whole broker. The relevant comparison includes
startup calibration, controller wakeups and calculations, measurement adapters,
pool resizing, and any lost throughput while the opening split converges. The
local `thread_broker_policy_overhead.py` therefore compares one release
binary in counterbalanced blocks at the same final split: oracle-pinned/no
broker, default responsive, sparse responsive, and freeze-after-convergence. It rejects a sample unless the
canonical RAD digest, counts, effective budget, final split, policy, convergence,
and CPU-accounting telemetry agree.

For every duration cell it reports paired absolute seconds and fractional ratios
for mapping wall, process wall, and aggregate process CPU. Formal acceptance
requires at least 30 complete blocks and a one-sided 95% bootstrap upper bound
<=1% for all three metrics on stable workloads. This percentage is deliberately
only a broad regression backstop. Process CPU is summed over all threads: at 64
fully occupied threads, 1% could hide an average 0.64 core, which is far too much
to call broker administration cheap.

The primary administrative-cost gate therefore directly sums lifetime
controller-thread and sampler-thread CPU. Each auxiliary thread reads its own
CPU clock once at entry and once at exit; there are no added clock reads in the
controller, sampler, decoder-task, or inflate loops. The runner reports this as
absolute nanoseconds, as a fraction of one core over mapping wall time, and as a
fraction of aggregate process CPU. Responsive and sparse-responsive modes must
have a one-sided 95% upper bound <=0.001 core (one millisecond of administrative
CPU per wall second). Freeze must use <=5 ms of administrative CPU through
convergence, after which both auxiliary threads are stopped and recurring cost
is zero. The aggregate <=1% comparison remains necessary to catch indirect
effects such as pool perturbation and CPU incurred on application threads.

Duration series must include approximately `{5 s, 1 min, 10 min}` real mapping
work. Freeze controller samples may grow by no more than 25% from shortest to
longest and its absolute CPU/wall delta must not grow materially; its fractional
overhead must decrease with duration. A responsive cadence is accepted only if
its chosen detection latency is explicit and both the direct and regression
backstop gates pass. Any policy must retain exact output and the same settled
allocation as the pinned control.

A real gate for a noisy, discrete controller decision must report its repetition
count, full range, final-allocation histogram, and outcome histogram alongside
the median. A single run is diagnostic only: it cannot qualify a mechanism that
can select different branches from noisy evidence. A cell that has exhibited
bimodality requires at least eight repetitions after a fix, with zero outcomes
outside its accepted allocation region and a per-run bound on startup-bracket
wall time.

The measured policy costs below make the distinction between broker
administration and allocation quality explicit. These are crossover-balanced
scATAC measurements on a 2-million-record, roughly 14--20 second workload; they
are evidence for policy selection, not a substitute for the still-optional
5-second/1-minute/10-minute duration matrix.

| policy | mapping wall versus pin-5 oracle | direct controller + sampler CPU | fraction of one core | important interpretation |
| --- | ---: | ---: | ---: | --- |
| pinned/no broker | 1.000 | none | none | measurement baseline |
| responsive, 25 ms steady probe | 1.217 | 4.965 ms | 0.02936% | administrative work is small; this older controller selected 6 rather than the pin-5 oracle |
| sparse-responsive, 5 s steady probe | 1.208 | 1.408 ms | 0.00835% | same selected split as responsive with 72% less administrative CPU |
| model-only freeze | 1.422 | 0.401 ms | 0.00197% | cheapest administration, but the one-point model selected producer 1 and is not suitable for this scATAC surface |
| freeze after full calibration | 1.016 median, 1.037 upper-95 | <=2.582 ms | <=0.0189% over 13.7 s | opening calibration found the oracle region before teardown |

A later eight-pair same-controller comparison isolates steady probing: sparse
versus 25 ms responsive had mapping ratio `0.9921` (upper-95 `1.0210`) and CPU
ratio `0.9941` (upper-95 `1.0196`), while median administrative CPU fell from
`4.902 ms` to `2.623 ms`. Thus sparse-responsive is the closest implementation
of the original low-overhead design: it preserves startup calibration and
recovery logic, then makes stable-run observation infrequent. Freeze is nearer
zero recurring cost, but only full-calibration freeze is safe on a non-monotone
response surface.

---

## 12. What would falsify this

- At a fixed pinned split, producer busy time misses Gate E's ground-truth error
  bounds, or accounted producer components explain less than 95% of CPU.
- Allocation-dependent costs vary by more than 10% and the observed fixed-point
  controller predicts held-out split throughput more than 10% incorrectly.
- A mapping-bound workload where the solved share is near zero but the measured
  optimum is materially above the producer floor — that would mean the pipeline
  has a cost the model does not see.
- A workload where `auto` lands >10% off a correctly-budgeted sweep. The claim is
  3–5%; twice that would put it in the range where the fixed splits it replaced
  were competitive.
- Convergence taking more than ~5 s or more than ~5 moves on any modality. DS2's
  three-step property is the thing being relied on.

---

## 13. Still open

Validation is deliberately tiered; “not in the default tier” is not the same as
“forgotten”:

| tier | default evidence | intended trigger |
| --- | --- | --- |
| light | thread-broker/unit checks, 32 position-balanced bulk-SE/Flex oracle runs, 16 short policy runs | routine controller feedback |
| normal | all five modalities at t8/t32/t64 (120 runs), 84-run 5 s/1 min/descriptive-10 min policy series, selected Gate G/TB-08 replay | default pre-merge qualification |
| comprehensive | 540-run five-pin oracle, 360-run policy matrix, full Gate E/G, short-duration grid, FIFO/mixed input, scATAC geometry/publication, cold-cache, and direct rapidgzip A/B | release/nightly or a change touching the associated risk |

The complete runnable catalog and exact local commands/fixtures live in
`completion-ledger.md`; those generated assets are intentionally ignored
by git. Gate E, the full real Gate G, FIFO/mixed input, scATAC geometry, and the
rapidgzip hot-path A/B have already passed. Comprehensive reruns remain opt-in
risk triggers, particularly measurement changes (Gate E and rapidgzip A/B),
lifecycle/input changes (FIFO), flow/cap changes (Gate G), and scATAC hot-path
changes (geometry/publication).

Two product questions remain outside normal gated completion:

- a true cold-page-cache end-to-end harness still needs a portable mechanism;
- two unreleased git dependencies block a crates.io release. The final
  rapidgzip changes are merged upstream at
  `276a41f77fb927e24cb0898a638a08b21eb048c6`; paraseq is still a feature branch.

Replacing the hand-written builders with `bon` remains optional maintenance,
not a release or correctness gate.

---

## References

- Kalavri, Liagouris, Hoffmann, Dimitrova, Forshaw, Roscoe. **Three steps is all
  you need: fast, accurate, automatic scaling decisions for distributed streaming
  dataflows.** OSDI'18. — true vs observed rate; the closed-form assignment; §2's
  catalogue of threshold-controller failures.
- Suleman, Qureshi, Khubaib, Patt. **Feedback-directed pipeline parallelism.**
  PACT'10. — `τᵢ = coresᵢ / avg_exec_timeᵢ`; LIMITER selection;
  revert-on-regression with a tried-set (§4.7).
- **FLIP-271: Autoscaling** (Apache Flink). — the production shape of a
  model-based autoscaler; `observed-scalability`; `detectIneffectiveScaleUp` and
  its default-off setting (§4.4).
- **FLINK-31215.** — non-scalable upstream bottlenecks defeating a rate model;
  the reason for §4.6.
- Kuo, Lim, Meerkov. **Bottlenecks in serial production lines: a system-theoretic
  approach.** Mathematical Problems in Engineering, 1996. — the signed indicator
  `e(d) = block₁(d) − starve₂(d)`, with an isolated zero where `starve` alone has
  a half-line. The formal statement of §2.
- **.NET CLR ThreadPool hill climbing** (`hillclimbing.cpp`; Hellerstein et al.'s
  Fourier band-pass design). — the in-flight sampling-bias guard behind §4.9.
- Page. **Continuous inspection schemes.** Biometrika, 1954. — CUSUM, implemented
  and then rejected in §4.8.
- Horstein. **Sequential transmission using noiseless feedback.** IEEE Trans. Inf.
  Theory, 1963; Waeber, Frazier, Henderson. **Bisection search with noisy
  responses.** SIAM J. Control Optim., 2013. — probabilistic bisection, held in
  reserve (§3.5).
- Combes, Proutière. **Unimodal bandits: regret lower bounds and optimal
  algorithms.** ICML'14. — OSUB, same.
