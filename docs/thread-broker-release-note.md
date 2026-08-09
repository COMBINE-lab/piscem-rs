# Thread-budget behavior change

Mapping commands now treat `-t` as one execution-slot budget shared by mapping
and parallel gzip decoding. The effective budget is capped by the process's
available parallelism and is recorded, with the requested value, in
`map_info.json`.

With `--decoder auto`, seekable gzip input uses a shared decoder pool and a live
thread broker adapts the aggregate split. Serial input and builds without the
`rapidgzip` feature give the full effective budget to mapping. `-t 1` always
uses the serial path.

Fixed controls have precise, changed semantics:

- `--decoder parallel=N` fixes `N` decode slots per decoder-capable input;
- `PISCEM_DECODE_SLOTS=N` fixes `N` aggregate decode slots and is the preferred
  control for reproducible oracle sweeps;
- legacy `PISCEM_RAPIDGZIP_THREADS=N` remains per decoder-capable input.

Any fixed control disables adaptation. Non-gzip files do not multiply a
per-file request. Oversized requests are visibly clamped to preserve one mapping
slot, zero/invalid/conflicting environment controls are rejected, and the
requested and applied allocations are emitted in `map_info.json`.

The broker reports final phase, requested and actual occupancy, empirical cap
reason, uncertainty, bounded trajectories, rejection evidence, and consumer
measurement components. Configuration errors still fail before mapping starts.
Runtime tuning failures (resize refusal/timeout, measurement startup, or broker
panic) are advisory: they retain structured failure telemetry, leave the last
valid under-budget split in force, and cannot suppress RAD finalization or mask
the mapping result.

Responsive mode offers an opt-in, self-bracketing opening policy for workloads
such as scATAC where stage cost depends on the allocation. It runs only when the
first stable model answer differs from the opening, measures at most three
200 ms candidate points over at most four seconds, and reports its point,
sample, wall, and outcome costs. Adaptive scATAC uses a 64-record
completed-progress publisher
during this decision phase; serial and pinned runs keep the generic 256-record
cadence. Busy-time clock
reads remain at 256 records even during calibration. Progress is written to
cache-padded, single-writer processor shards with relaxed cumulative stores, so
the finer cadence passed the formal <=1% wall- and CPU-overhead gates.
Stable scATAC runs probe every 5 s after convergence by default; a builder or
the validation environment hook can restore the 25 ms responsive cadence.
Freeze-after-convergence is the explicit model-only low-overhead policy: it
skips opening calibration and stops recurring controller, sampler, and
fine-publication work once settled. Freeze-after-full-calibration first runs the
bounded opening bracket, then performs the same teardown; it is the
appropriate freeze policy where the one-point model can miss the response
curve.

The model answer remains primary during the opening bracket. An apparent loss
at the short opening horizon is extended to the ordinary ratification sample
count before exploration; the common accept path still stops at the short
horizon. A regressed or inconclusive comparison may trigger an adjacent response
point, but that point must beat the higher measured opening/model rate. Only a
greater-than-5% gain with non-overlapping intervals can replace the model. The
opening is a pivot and high-water evidence bar, never the fallback, and
deadline-limited confirmation likewise retains the model.
Independent persistent cap evidence can confirm the model target without local
exploration: if the producer demonstrably cannot use slots above that target,
short-term throughput disagreement is not evidence for buying them back.

Measured on a crossover-balanced 14--20 second scATAC workload, direct
controller-plus-sampler CPU was 4.965 ms for 25 ms responsive monitoring,
1.408 ms for 5 s sparse-responsive monitoring, and 0.401 ms for model-only
freeze. A later same-controller eight-pair comparison reduced median
administrative CPU from 4.902 to 2.623 ms with sparse monitoring, while its
mapping ratio versus 25 ms monitoring was 0.992 (upper-95 1.021). Model-only
freeze was cheapest but chose a poor split on this non-monotone surface;
full-calibration freeze instead stayed within 1.016 median/1.037 upper-95 of
the oracle and used at most 2.582 ms of administrative CPU. Sparse-responsive
is therefore the closest policy to the original low-overhead responsive design.
Adaptive scATAC now opens with four aggregate decoder slots at every budget and
leaves the true shared-pool safety floor at one. Four is only a hint: when the
model disagrees, the bracket measures its answer and at most two points adjacent
to the opening. This lets the same mechanism select producer five on the
negative-scaling t8 surface and producer one on the mapping-heavy t32/t64
surface, with no `budget <= 8` policy threshold. Serial and pinned controls are
unchanged. scRNA and Flex retain their measured quarter-budget opening and the
default fixed opening policy, so they perform no bracket work.

The oracle runner also refuses to label a boundary winner as an oracle. A best
fixed point must have measured allocations on both sides; otherwise the cell is
reported as an unbracketed candidate and fails qualification.

On the final 2-million-read scATAC fixture, two t8 runs selected producer five
in two points and 16 samples. Their 14.248 s median was 1.036 times the 13.7575 s
pin-five oracle median; bracket wall was 1.353--1.453 s and controller CPU was
0.404--0.465 ms (0.00305% of one core at the median wall). At t32 the bracket
selected producer one in one point/two samples: 7.938 s mapping, 0.800 s bracket
wall, and 0.205 ms controller CPU. At t64 two adaptive runs likewise selected
one in one point; their 4.380 s median was below the 4.925 s same-binary pin-one
median, although both very short arms were noisy. Every output had the same
canonical digest. These are startup-bracket measurements; the settled-policy
overhead measurements above remain the evidence for recurring cost.

A follow-up review found that the single t32 run hid a bimodal 5/8 model versus
3/8 opening result. The corrected final binary (`78cebc9c...`) makes the opening
a candidate-placement pivot rather than a fallback, requires alternatives to
beat the higher opening/model evidence bar, and accepts an independently
confirmed empirical cap. Eight t32 repetitions then selected producer one in
8/8 runs: 4.526--11.859 s, median 7.808 s (0.907 times the 8.61 s pin-one
reference), with seven `model_selected` outcomes and one deadline-limited model
fallback. Bracket wall was 0.700--1.001 s normally and 4.000111 s at the single
deadline (111 microseconds of wakeup overshoot); controller CPU was
0.129--0.363 ms. Two t8 guards selected producer five in 2/2 runs, with
1.052--1.153 s bracket wall. All ten adaptive outputs and two fixed controls
retained the canonical digest above.

The exact cumulative rapidgzip signal is feature-gated upstream; its disabled
build compiles the hot-path accounting out, and its enabled build passed the
paired no-measurable-overhead gate. Release remains blocked on the unreleased
dependency pins documented in the completion ledger.
