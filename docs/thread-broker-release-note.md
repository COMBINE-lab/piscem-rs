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
measurement components. Resize, timeout, and sampler failures now fail the
mapping command instead of returning a successful zero/default report.

Responsive mode offers an opt-in bounded nonlinear response fallback for
workloads such as scATAC where adding mapping workers can reduce throughput;
piscem enables it only for scATAC. scATAC uses a
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
The adaptive scATAC opening is the midpoint rather than the generic
mapping-biased split, providing a safe anchor for that response search without
changing serial or pinned controls.

Release remains blocked until the cumulative producer-measurement, full
modality/oracle, overhead, and upstream-dependency gates in
`THREAD-BROKER-AUDIT.md` pass.
