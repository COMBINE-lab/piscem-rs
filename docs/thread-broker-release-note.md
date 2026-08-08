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

Responsive mode offers an opt-in bounded nonlinear response fallback for
workloads such as scATAC where adding mapping workers can reduce throughput.
Piscem enables that search for small-budget scATAC and exposes an application
cap on the explored producer allocation. scATAC uses a
64-record completed-progress publisher only during this decision phase and
returns to the generic 256-record cadence after convergence. Busy-time clock
reads remain at 256 records even during calibration. Progress is written to
cache-padded, single-writer processor shards with relaxed cumulative stores, so
the finer cadence passed the formal <=1% wall- and CPU-overhead gates.
Stable scATAC runs probe every 5 s after convergence by default; a builder or
the validation environment hook can restore the 25 ms responsive cadence.
Freeze-after-convergence is the explicit model-only low-overhead policy: it
skips nonlinear probing and stops recurring controller, sampler, and
fine-publication work once settled. Freeze-after-full-calibration first runs the
bounded nonlinear/local search, then performs the same teardown; it is the
appropriate freeze policy where the one-point model can miss the response
curve.

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
Adaptive scATAC opens with four aggregate decoder slots, the stable region in
the measured fixed-split surface. At larger budgets the ordinary cost model can
still grow above four, but the optional nonlinear startup search is disabled
because its noisy probes did not amortize. Serial and pinned controls are
unchanged. scRNA and Flex use a measured quarter-budget opening.

The exact cumulative rapidgzip signal is feature-gated upstream; its disabled
build compiles the hot-path accounting out, and its enabled build passed the
paired no-measurable-overhead gate. Release remains blocked on the unreleased
dependency pins documented in the completion ledger.
