# Threading, decompression, and the tiny-dict prefilter

Design notes for the mapping pipeline's use of threads and its k-mer prefilter.
Everything here is stated with the measurement that justifies it; where a
plausible-sounding alternative was tried and lost, that is recorded too, so it
is not retried from first principles.

Reference workloads used throughout:

| name | index | k-mers | reads | gzip structure |
|---|---|---|---|---|
| **Flex** | 10x probe panel | 1.5 M | 150 bp | multi-member, ~186 KB |
| **PBMC** | gencode v49 | 96.3 M | 91 bp | multi-member, ~12 KB |
| **bulk** | gencode v49 | 96.3 M | 150 bp | single-member |

The bulk set was ERR3239276, which **maps at only 11%** against gencode v49
(measured at two different offsets in the file, so not a head artifact). It is
fine for I/O and threading work, where only the byte volume matters, but it is
the wrong choice for anything that depends on lookup outcomes. Use
`SRR21186103` instead: 150 bp, **88.8%** mapping rate on the same index.

`SRR21186103` is available as a **pair**, and both mates are needed to cover the
paired-end path (`map_pe_fragment`, two mates filled per group under separate
per-file mutexes). Both are at
`/scratch1/rob/rshash_testing/human_reads/SRR21186103_{1,2}.fastq.gz`, fetched
from ENA and byte-identical to the archive copies; both are gzip-family
(`MarkerWindow`, zero `Z_SYNC_FLUSH` markers). Paired, the run maps 70.3% of
26.1 M pairs.

The pair earns its place by being the one reference workload that is
**mapping-bound**: 10.0% consumer idle at `-t 32`, where Flex measures 93.8% and
bulk SE 50.7%. Any thread-allocation policy validated only against decode-bound
inputs can pass by always favouring decode, so a workload where the right answer
is *serial* is what keeps that honest.

Machine: AMD EPYC 9575F, 2 sockets x 64 cores.

## 0. `pigz` output and rapidgzip (fixed upstream in 0.2.1)

**Resolved — kept because the failure was invisible and is easy to re-encounter
on an older dependency.** Against `rapidgzip-core` **0.2.0**, `pigz`-written gzip
was never admitted to the parallel marker path: it committed to
`DecoderPath::Sequential` during admission and decoded single-threaded for the
whole file, whatever the thread budget or the calibrator's verdict.
`gzip`- and `zlib`-written streams of identical content were admitted.

Reported as
[COMBINE-lab/rapidgzip-rust#29](https://github.com/COMBINE-lab/rapidgzip-rust/issues/29)
and fixed in **0.2.1**, which this crate now requires. Every fixture is
`MarkerWindow` or `DenseMembers` again; measured on the same files, `bulk8M`
went from ~1770 MB/s sequential to 6493 MB/s.

The cause, for anyone meeting a similar symptom: at a `pigz` flush point the
exact decoder reports the boundary *before* a short zero-output block while the
speculative searcher resolves to the boundary *after* it, so
`run_admission_wave`'s strict `chunk.start_bit != previous_end` equality
rejected. The cleanest instance was a 35-bit delta — a 3-bit block header plus a
32-bit `LEN`/`NLEN`, i.e. an empty stored block (`Z_SYNC_FLUSH`). Neither
single-member-ness nor byte alignment was the discriminator; admitted `gzip`
boundaries are also non-byte-aligned and resolve exactly.

**What survives the fix.** A correct calibrator can still select `parallel` and
have the decoder report `Sequential`, because that is a property of the input's
encoding rather than a bug. `io::decode_budget` logs exactly that mismatch at
`WARN`, which is the only way a user learns their file is the limiter — the
decoder otherwise reports `Idle` with `spawned_workers = 0`, indistinguishable
from having nothing to do. Diagnosing a suspected case takes one line:

```bash
head -c 50000000 reads.fq.gz | python3 -c \
  "import sys; print(sys.stdin.buffer.read().count(b'\x00\x00\xff\xff'))"
```

Nonzero means `Z_SYNC_FLUSH` markers, i.e. `pigz`-family output. On 0.2.1 that is
no longer a problem; on anything older it is decisive.

## 1. Decompression is usually the constraint, not mapping

paraseq holds one mutex per input file, so exactly one thread can inflate a
given file at a time — a hard ceiling of two inflate streams for a paired run,
regardless of `-t`. On the full Flex archive (2 x 149.5 GB, 2.30 B read pairs)
that capped the run at ~11 min with mapping threads ~80% idle, and doubling
`-t` from 16 to 32 changed wall clock by 1%.

`rapidgzip-core` (feature `rapidgzip`) decodes a *single* gzip member in
parallel, lifting the ceiling without changing paraseq's structure:
**11.28 min -> 1.86 min**, output bit-identical.

Non-gzip input (plain, zstd, bz2) keeps the niffler path and allocates no
decoder threads.

## 2. Current decoder selection and allocation

`-t` is the total execution-slot budget shared by mapping and parallel gzip
decoding. The effective budget is capped by `available_parallelism`; the
requested and effective values are both written to `map_info.json`. Serial
input, builds without `rapidgzip`, and an effective budget of one give the whole
budget to mapping.

`--decoder` has four forms:

- `auto` opens seekable gzip through the shared decoder pool and lets the live
  thread broker adapt the aggregate mapping/decode split;
- `serial` disables the parallel pool and gives mapping the full budget;
- `parallel` forces the parallel path where the input supports positional reads,
  but still lets the broker adapt the split;
- `parallel=N` fixes `N` decode slots per decoder-capable input, reconciles that
  request after the files are opened, and disables adaptation.

Every fixed request is clamped visibly when necessary to preserve at least one
mapping slot. `PISCEM_DECODE_SLOTS=N` is the aggregate fixed control intended
for oracle sweeps. The legacy `PISCEM_RAPIDGZIP_THREADS=N` remains a per-file
fixed control. Setting both environment controls is an error. Non-gzip paths do
not count as decoder-capable inputs, and FIFOs/process substitutions always use
the serial path because they cannot be rewound or positionally read.

In adaptive mode, one shared `ExecutionPlan` sizes both pool maxima, the opening
split, and the broker budget. The broker measures non-blocked mapping and decode
work during the real run, rejects windows whose inter-stage buffer is materially
filling or draining, and applies a solved split using shrink acknowledgement
before growth. Resize failures and sampler failures abort with typed errors.
Reports include the requested and actual split, final phase, uncertainty,
empirical cap reason, bounded trajectories, rejections, controlled/auxiliary
occupancy, and consumer component timings in `map_info.json`.

Effective budgets up to eight use a 25 ms cadence and 100 ms warm-up so short
jobs can adapt after the observed startup transient. They open with two decode
slots where the budget permits. A producer that still has runnable work queued
may veto a model-requested shrink, but pressure never grows or sizes the split;
this is the bounded fallback for allocation-dependent decoder scaling.
scATAC instead uses one four-slot opening hint at every budget and opts into a
bounded startup bracket. Its safety floor remains one: if the first stable model
answer differs from four, the broker measures that answer and at most two local
alternatives rather than encoding a measured optimum as an unreachable floor.
Serial and fixed allocations retain their exact requested split.

The decoder busy-time signal is integrated from lock-free per-decoder executing
worker counts. While a scheduling decision is open, a dedicated sampler uses a
jittered cadence averaging 3 ms and trapezoidal integration; after the broker
settles it drops to 25 ms and returns to the denser cadence only for a resurvey.
Sample counts and observation time are reported in `map_info.json`. Generated
burst/cadence accuracy tests and same-binary sampler-on/off overhead controls
guard this design.

Applications can select steady-state behavior with the thread-broker builder.
`SteadyStatePolicy::Responsive` is the default and preserves regime-change
adaptation. Applications may opt into `OpeningPolicy::Bracket` when
allocation-dependent stage scaling can invalidate the one-point service-cost
model; piscem does so for scATAC. The startup-only experiment is triggered by
model/opening disagreement, is bounded by point count and wall time, and reports
its cost and outcome. Other modalities keep `OpeningPolicy::Fixed` and pay no
bracket cost. Its `steady_probe_interval` is
independently configurable: it does
not weaken startup calibration or ratification, and the sampled decoder adapter
scales toward roughly four observations per probe without exceeding its measured
25 ms accuracy floor. Stable scATAC uses a 5 s responsive interval by default;
`PISCEM_THREAD_BROKER_PROBE_INTERVAL_MS=25` restores normal-resolution
monitoring without changing startup calibration. Other modalities retain their
ordinary cadence. `SteadyStatePolicy::FreezeAfterConvergence` is the cheaper
model-only path: it skips opening calibration, then terminates the
controller and producer sampler after guarded model convergence. The separate
`SteadyStatePolicy::FreezeAfterFullCalibration` runs responsive mode's bounded
opening calibration before the same teardown. Both freeze modes
have no recurring broker work and cannot react to a later workload change; only
the full-calibration variant is suitable when the one-point model is known to
miss the useful response curve.

On the long scATAC validation fixture, an eight-pair direct comparison found no
throughput penalty from the sparse default: one-sided 95% upper sparse/25 ms
ratios were 1.0211 for mapping wall, 1.0209 for process wall, and 1.0197 for
aggregate CPU. Median controller wakeups fell from 566.5 to 144.5 and decoder
monitoring observations from 431 to 10.

scATAC publishes completed-record mapping progress every 64 records only while a
responsive decision is open; its more expensive records made the generic
256-record cadence too coarse for 25 ms windows. Collecting more coarse windows
does not repair that quantization, so the implementation increases publication
resolution during calibration instead. It returns to 256 in responsive steady
state, and freeze reduces new batches to drop-only publication after convergence.
Only the item counter uses the 64-record cadence; busy-time clock reads remain at
256 records. Each mapping processor owns a cache-padded progress shard and is
the shard's only writer, so publication is one relaxed cumulative store rather
than a contended shared increment. The controller sums the shards only at its
sampling cadence. Generic consumers retain the shared 256-record meter.

The validation-only `PISCEM_SCATAC_PROGRESS_FLUSH_EVERY=N` hook permits
same-binary 64-versus-256 fixed-split tests. The formal 30-pair,
2-million-record comparison was exactly balanced (15 fine-first and 15
generic-first), retained canonical output and counts in all 60 runs, and passed
the one-sided 95% <=1% overhead gates. Position-stratified upper ratios were
`1.00419` for mapping wall, `1.00434` for process wall, and `1.00441` for
aggregate process CPU; paired medians were `0.99933`, `0.99895`, and `0.99965`.
The local measurement harness and raw-evidence location are inventoried in the
ignored completion ledger.

Piscem's same-binary A/B hooks are
`PISCEM_THREAD_BROKER_POLICY=responsive|freeze-after-convergence|freeze-after-full-calibration` and
`PISCEM_THREAD_BROKER_PROBE_INTERVAL_MS=N`. Machine telemetry records the
selected policy, effective probe interval, controller samples/lifetime, and
whether monitoring stopped at convergence. These hooks are intended for policy
validation; fixed oracle runs still use `PISCEM_DECODE_SLOTS=N` and disable the
broker. The policy-overhead runner compares fixed, default responsive, sparse
responsive, model-only freeze, and full-calibration-freeze modes selected by its
manifest, using paired absolute and fractional CPU/wall cost including
minute-scale runs. Its <=1% process-CPU comparison is an aggregate
regression backstop, not the administrative budget: on 64 saturated threads 1%
would allow 0.64 core. Controller and sampler threads therefore also record
their complete lifetime CPU using one clock read at entry and one at exit. The
direct responsive gate is <=0.001 core (1 ms CPU per wall second); freeze is
limited to <=5 ms through convergence and then stops both threads.

An optional rapidgzip feature also records worker and auxiliary *thread-lifetime
CPU time* for validation. It reads the CPU clock only when a thread registers
and exits and updates a cumulative counter only at exit; it adds nothing to a
decode-task or inflate hot loop. The feature is disabled by default and is not
the controller signal. It passed separate direct-decoder and application-level
feature-on/off gates whose one-sided 95% wall- and CPU-overhead bounds were at
most 1%. See `THREAD-BROKER-AUDIT.md` for the remaining real-workload gates.

<details>
<summary>Archived pre-run probe and supervisor experiments (superseded)</summary>

The material below records experiments that motivated the live broker. It is
historical and does not describe the current command-line behavior.

### Choosing the decoder: threads per file, not file count

The serial (niffler/libz-rs) path opens exactly one inflate stream per input
file, so its supply is `files x per-stream rate` and does **not** respond to
`-t`. Demand is `map_threads x per-thread consumption`. The parallel decoder
therefore pays off on the *ratio*, and a file-count threshold is only right at
the thread count it was calibrated at.

Measured on the 10x Flex probe panel (multi-member archives), wall-clock
speedup of the parallel decoder over the serial one:

| files | threads/file | speedup |
|------:|-------------:|--------:|
| 2 | 4 | 1.09x |
| 2 | 8 | 1.92x |
| 2 | 16 | 3.05x |
| 4 | 4 | 1.07x |
| 4 | 8 | 1.54x |
| 4 | 16 | 2.08x |
| 8 | 4 | 0.92x |
| 8 | 8 | 1.20x |

The serial path's wall clock floors at a level set purely by file count -- ~9.8 s
at one pair, ~5.0 s at two, ~3.0 s at four, regardless of `-t` -- and the
parallel decoder wins exactly where threads push demand past that floor.

The previous rule (`files >= 8` -> serial) switched on at 8 files / 32 threads
where the parallel decoder *loses* (0.92x), and off at 8 files / 64 threads
where it wins (1.20x). Replaced by `map_threads / files >= 4`. Verified:

| pairs | -t | ratio | serial | old rule | new rule |
|---|---|---|---|---|---|
| 2 | 32 | 8 | 5.05 | 3.22 | 3.22 |
| 2 | 64 | 16 | 5.21 | 2.53 | 2.51 |
| 4 | 16 | 2 | 5.23 | 5.23 | 5.27 |
| 4 | 64 | 8 | 2.95 | 2.93 | **2.55** |

**Both rules above are historical.** Neither a file-count threshold nor a
threads-per-file ratio survives, because both stand in for the term that actually
decides: per-thread mapping cost, measured at 0.064 GB/s on a 96.3 M-k-mer
transcriptome against 0.43 GB/s on a 1.5 M-k-mer probe panel. A 6.7x spread means
the *correct* ratio is about 3.5 threads/file for cheap mapping and about 37 for
expensive, and no constant spans that. The table is kept because it is the
crossover surface any replacement has to be checked against.

The payoff asymmetry does survive, and still sets which way ties break: enabling
the parallel decoder needlessly measured -4.9% wall and +2.1% CPU on a
transcriptome where decode never binds, while failing to enable it when needed
costs up to 3x wall. The cost of the wrong "yes" is threads; the cost of the
wrong "no" is time.

### Calibration: run the pipeline, do not model it

`io::probe` decides by **running the real pipeline briefly** — the same
`Collection` API, reader type, batch size, per-file mutexes and thread
distribution production uses — and measuring what fraction of mapping-thread
capacity is not spent mapping. It is deliberately independent of the mapping
kernel (it takes a `ProbeKernel`) so downstream consumers such as salmon can
share it.

**An earlier design compared two separately-measured rates and is gone.** Both
rates were accurate — the producer estimate came within 2.6% of a whole-file
measurement of identical work, the consumer within 3% — but validated against
the crossover surface above, `decide` got **4 of 8 points wrong**, every one
predicting serial where the parallel decoder measurably won, by up to 1.92x. The
flaw was structural: it modelled supply as `files x rate-measured-alone`, but the
thread holding the per-file mutex also maps, so achieved supply falls exactly as
`-t` rises — the model was most optimistic where it was most wrong.

Two traps found while building the replacement, both worth remembering because
each produced confident, plausible, inverted numbers:

- **Per-batch time accounting cannot be sampled.** `READER_BATCH_SIZE` is 16384,
  so at a realistic 5 us/record a batch takes ~80 ms against a 25 ms window.
  Windows containing no completed batch read a delta of zero, and the probe
  reported **100% starved for a workload mapping flat out**. Mapping time is now
  published every 256 groups.
- **The verdict landed on the startup transient.** Every workload stopped at
  exactly `min_windows` reporting ~99% starved, because at 75 ms `paraseq` is
  still spawning workers and nothing has been mapped yet — so *every* input looks
  decode-bound. A 150 ms warm-up is now discarded and re-based.

With both fixed the measure is monotone across a 16,000x range of synthetic
consumer cost, from 99.6% starved to 1.1%.

**Forcings** still apply first, and only where the answer cannot depend on
measurement: no regular file among the inputs, or a single mapping thread. The
threads-per-file bound that used to live here is gone — it encoded consumer cost
as a constant, which was empirical rather than logical, and consumer cost is
exactly what the probe measures.

### Never probe what cannot be rewound

The probe reads a prefix and the run re-opens from the start. On a FIFO those
bytes are **gone** — there is no second read — so probing one silently drops
reads from the output. That is data loss, not a slowdown, and it is the one
failure mode here that corrupts results rather than costing time.

The guard is per *group*, because `paraseq`'s paired `fill` reads both mates
together: a regular R1 beside a FIFO R2 cannot be probed without consuming the
pipe. Groups containing anything non-regular are excluded; if none remain the
probe is skipped entirely and the run keeps its default.

The downgrade is **per file**, not per run: a FIFO among the inputs does not cost
the regular files their parallel decoder, since `open_input` re-checks each path.
Under `--decoder auto` the partial downgrade is reported at `INFO`; under an
explicit `parallel` it is a `WARN` naming the offending paths, because "parallel
was requested but you got serial" is actively misleading when one file of eight
was affected.

Tested by construction: opening a FIFO with no writer blocks forever, so the
tests create writer-less FIFOs and a completing test *is* the proof that nothing
opened them.

### Overriding the choice

`--decoder <auto|serial|parallel|parallel=N>`, defaulting to `auto`. For someone
who knows something the probe cannot: a slow network filesystem, a shared node
where spending cores on decode is antisocial, or reproducing a measurement.
`parallel=N` sets the per-file worker count as both ceiling and starting point —
naming a number means wanting it used, not ratcheted toward.

A preference outranks measurement but **not** the forcings: `parallel` on a FIFO
still yields, because the parallel decoder degrades to sequential there and a
preference cannot make an input seekable.

**A request that cannot be honoured warns**, at `WARN`, visible without setting
`RUST_LOG` — the subscriber floors at `warn` when the variable is unset, since a
message saying a flag was overridden is useless if only visible to someone who
already suspected a problem. Two cases: `parallel` on a non-regular input, and
any explicit decoder on input that is not gzip, where neither path competes.

### scATAC

`map-scatac` opened every input through the niffler-only helper, so `-t` and
`--decoder` did nothing for it whatever their value. It now goes through the
same `plan_thread_budget` / `open_input` path as the other two.

Measured, and it is a consistency fix rather than a speedup. 10x `atac_pbmc_5k`
lane L001 (107.8 M reads, 3 files) against a k=23 chr1+chr2 GRCh38 index, `-t 32`,
3 reps each:

| | wall | CPU |
|---|---|---|
| serial | 324.2 s | 10,310 s |
| parallel | 327.3 s | 10,426 s |

+1.0% wall and +1.1% CPU for the parallel decoder — inside the noise floor and
if anything marginally worse. scATAC opens **three files per sample**, so the
serial path already supplies three inflate streams, and mapping a 397 M-vertex
index is costly enough per read that decode never binds: the supervisor reported
`peak busy 3` against a budget of 16, having found no starvation to answer.

The caveat runs against the change rather than for it: a two-chromosome
reference makes mapping cheaper than a whole-genome run and so biases the
workload *toward* decode-bound. Even so the parallel decoder does not win, so on
full GRCh38 the case would be weaker still. The value here is that the policy
now applies at all and reaches the right answer, not that the answer is faster.

### Non-regular inputs (FIFOs, process substitution)

`rapidgzip` gates its parallel path on `file_type().is_file()` and falls back to
*sequential* decoding otherwise, so it has nothing to offer on a pipe. Such
inputs take the serial path unconditionally (`calibrate::Reason::NonSeekableInput`).

This also fixes a hang. `open_gz_rapidgzip` sniffs the gzip magic by opening the
path, reading two bytes, and closing it, then calls `decoder.open(path)` which
opens it *again*. On a FIFO the sniff consumes those bytes and closes the read
end, and the re-open blocks forever waiting for a writer that has exited — so
`-r <(zcat ...)` hung indefinitely with the `rapidgzip` feature enabled, while
the serial build handled it fine.

Confirmed to be ours and not upstream: `rapidgzip-core` alone, opening the same
FIFO with no prior open of the path, decompresses it correctly and byte-for-byte
identically to the regular file. Nothing was reported upstream.

## 3. Thread budget: how `-t` is split

**`-t` is now the whole budget**, divided into mapping threads and a decode
budget. It previously meant mapping threads *alone*, with decode allocated on
top, so a run used about 1.5x the cores it named — which is what cgroup-limited
runs tripped over. Existing timings shift accordingly.

**Allocation matters far more than budget size.** At a fixed 64 threads the split
was worth **5.75x** (11.04 min with every thread mapping and no parallel decode,
versus 1.92 min balanced), while doubling 64 -> 128 bought only 9%. It was
nevertheless a constant — `decode = map / 2` — for as long as this module
existed, which made the one parameter that mattered the one thing never
measured.

A constant cannot work, because the right answer moves an order of magnitude with
the index: the knee sat at 16 workers/file for Flex, 4 for bulk, and effectively
0 for PBMC on a 96.3 M-k-mer transcriptome. So `plan_thread_budget` now derives
it from two measured rates: per-thread mapping consumption, and achieved decode
supply. Throughput is `min(demand, supply)`, so it evaluates every split and
keeps the best.

### A decoder worker is not worth a stream

The obvious model — a worker supplies what one serial stream does — is wrong, and
measurably so. Measured on gencode v49, one gzip file, `-t 32`, best of 3,
**mapping seconds only** (whole-run wall clock buries this under 1.58 s of index
load):

| decode threads | mapping |
|---:|---:|
| 1 (what the closed form chose) | 5.70 s |
| 2 | 5.68 s |
| 4 | 3.90 s |
| **6 (optimum)** | **3.57 s** |
| 8 (what the two-point split chooses) | 3.87 s |
| 12 | 4.40 s |

The closed form lands **60% off**, always toward too little decode. The reason is
Amdahl, not tuning: `paraseq` holds one mutex per file across `fill`, which
**decodes *and* parses**. Workers parallelise only the decode half; parsing stays
serialised per file whatever the worker count. A worker therefore buys a fraction
of a stream — about `0.28 * s` here — and the fraction depends on the file's
structure, since `DenseMembers` and `MarkerWindow` parallelise differently.

No constant fixes that, so the probe measures supply a **second time** with
decoders running, and the split interpolates between the two points and holds
flat past the second (supply saturates; extrapolating promises throughput no
worker count can deliver). That lands within **7%** of the optimum — above the
~1-3% noise floor, and the residual comes from linear interpolation across a
concave curve, which biases the choice high.

Cost: a second probe, a few hundred ms, paid only when the first probe finds
enough idleness to be worth measuring.

### Serial versus parallel is a throughput comparison, not a threshold

`parallel_gain` compares predicted throughput at the best split against spending
nothing on decode, and the parallel decoder is used when that clears a few
percent. A starvation threshold cannot answer this: idleness says the mapping
threads are waiting, not that decoder threads would fix it — a file rapidgzip
decodes sequentially is 100% starved and cannot be helped at all. It would also
have to be calibrated, which is the fitted constant this design removed.

`STARVATION_SCREEN` (0.05) survives only as a cheap gate on whether the second
probe is worth taking, and it is derived rather than fitted: idle fraction `f`
caps recoverable speedup at `1/(1-f)`, and enabling parallel decode needlessly
measured -4.9% wall / +2.1% CPU, so below ~5% idle nothing can pay.

Observed, `-t 32`:

| workload | files | idle | gain | decision |
|---|---:|---:|---:|---|
| Flex probe panel | 1 | 78.5% | 2.60x | parallel |
| gencode v49 | 1 | 51.3% | 1.88x | parallel |
| gencode v49 | 8 | 4.5% | — (screened out) | serial |
| scATAC, 10x pbmc_5k | 3 | 0.0% | — (screened out) | serial |

### Growth is one-way within a run

Decoders open at their ceiling and `io::decode_budget` only ever raises it.
Shrinking was tried and oscillates with a ~1 s period: at 4 workers the decoder
is `DecoderBound` so the share grows to 8; at 8 it gets ahead of the parser and
reports `ConsumerBound` so the share halves; forever. Upstream's own controller
already adapts downward — measured directly, a converged decoder drops from 16
active workers to 1 under `ConsumerBound` with no help from here — so shrinking
the ceiling is a second controller fighting the first.

The supervisor's remaining job is the one upstream cannot do: dividing a fixed
shared budget across the decoders a multi-file run opens. It accounts in
**live workers**, not granted ceilings. Summing ceilings deadlocked it outright
once decoders stopped being clamped at open: every decoder reported
`worker_limit == per_file_ceiling`, roughly twice the fair share, so `held` came
to about twice the budget, `spare` was zero on every sample, and the growth
branch never fired.

### scATAC

`map-scatac` now takes `--decoder` and calibrates like the other two, rather than
opening through a niffler-only helper where `-t` and `--decoder` did nothing.

It is a consistency fix, not a speedup, and the calibrator agrees: 10x
`atac_pbmc_5k` lane L001 against a k=23 chr1+chr2 index at `-t 32` probes
**0.0% idle** and chooses serial. scATAC opens three files per sample, so the
serial path already supplies three inflate streams, and a 397 M-vertex index
makes mapping costly enough per read that decode never binds. Measured
end-to-end the parallel decoder cost +1.0% wall and +1.1% CPU. The value is that
the policy now applies at all and reaches the right answer by measurement.

</details>

## 4. The tiny-dict prefilter

A blocked Bloom filter in front of the hashbrown probe, gated on two properties.
See `tiny-dict/src/prefilter.rs` and the `FilterGate` in `tiny-dict/src/lib.rs`.

**It is consulted only on genuine searches.** `try_extend_fw` (within-unitig
extension) and piscem's fast-positive skip check both resolve against the SPSS
and return before it. So the governing statistic is the *search*-miss rate,
which differs sharply from the per-k-mer rate because extensions absorb hit runs:
Flex 99.0% (vs 79.12% per k-mer), PBMC 86.6% (vs 45.50%).

**Gate A — cache residency (8 MB).** Not a memory bound. Benefit tracks index
size, not miss rate: -16.3% CPU on the 1.5 M-k-mer panel versus -5.1%/-3.8% on
the 96.3 M-k-mer transcriptome despite comparable search-miss rates.

**Gate B — observed search-miss rate.** Sample 200 K searches, enable at >= 70%.
Gate A is the real discriminator (PBMC's 86.6% would pass B); B covers the case
A cannot see — a small index queried by reads that mostly map.

### Do not "improve" the filter without measuring first

Three changes were made and **all reverted** after measurement:

1. Multiply-shift range reduction replacing power-of-two masking. Correct in
   itself (masking inflated an 8 bits/key request to 11.2 allocated, ~40%
   waste) — but the waste was silently buying FPR.
2. 512-bit cache-line blocks with 5 bits/entry, fixing FPR 7.39% -> 2.34%.
3. Parallel construction (3.118 s -> 0.126 s at 96.3 M keys), which then
   justified relaxing Gate A.

Result: **strictly worse on every workload** — flex +8.3% CPU (from -16.3%),
pbmc +16.6%, bulk +6.2%. Both mistakes were the same shape, optimising a
quantity that was not the constraint:

- **FPR is nearly irrelevant here.** At 99% search-miss a miss is rejected
  either way, and a false positive costs only the probe that would have
  happened anyway. What is paid on *every* search is probe cost, and the
  rewrite tripled it.
- **Gate A was never about construction time.** It is about cache residency.

### External crates are not faster here

Benchmarked against `fastbloom` 0.17 ("the fastest Bloom filter in Rust") on the
real access pattern — 1.5 M keys, 99% absent queries — with a fast hasher
(rapidhash) given to both:

| | memory | build | probe | FPR |
|---|---|---|---|---|
| ours | 2.1 MB | 0.003 s | **1.68 ns/op** | **1.388%** |
| fastbloom | 1.5 MB | 0.007 s | 11.46 ns/op | 2.159% |

Not a knock on the crate — it is generic (`impl Hash` dispatch, `BuildHasher`
boundary, runtime-configured sizing). Our path is a specialised leaf: keys are
already `u64`, the hash is one `mulx`, the probe is one load and a branchless
compare. Profiling puts the filter probe at 19.4% of cycles, so the ~10 ns/op
difference is worth roughly 2% of whole-run time.

Xor/binary-fuse filters are also the wrong trade: 3 probes where we need 1, and
peeling construction of ~10 s at 96 M keys, to improve an FPR that barely
matters.

**If the remaining 19.4% is worth attacking, the lever is fewer filter probes
(the skip machinery issuing fewer searches), not a faster probe.**

## 5. A negative prefilter does not pay on the sshash backend

Tried and **reverted**. The motivation looked strong: dictionary seeds are
overwhelmingly misses, because a hit is absorbed into a within-unitig extension
run while every miss costs its own seed.

| workload | mapping rate | seeds / 500K reads | seed miss rate |
|---|---|---|---|
| 10x Flex R1 | 97.2% | 46.3 M | 97.7% |
| ERR3239276 | 11.1% | 53.7 M | 96.9% |
| SRR21186103 | 88.8% | 10.0 M | 82.2% |

A blocked Bloom filter over the k-mer set, probed before the seed chain, was
built and swept over 1/2/4/8/16 bits per key. Output stayed byte-identical at
every size. It lost everywhere, badly:

| bits/key | filter size | SRR21186103 CPU | Flex R1 CPU |
|---|---|---|---|
| off | -- | **30.52 s** | **48.82 s** |
| 1 | 16 MB | -- | 77.73 s |
| 2 | 32 MB | 47.85 s | 87.59 s |
| 8 | 128 MB | 43.15 s | 78.65 s |
| 16 | 256 MB | 47.87 s | -- |

Construction was not the cost -- it parallelizes to 0.09-0.12 s. Three reasons
it fails, and they compound:

- **The bucket cache already removed the work it was meant to avoid.** 71% of
  seeds now resolve their bucket from a cached value. The filter adds a probe in
  front of a chain that is no longer expensive.
- **Hashing destroys locality.** Consecutive seeds within a read touch nearby
  unitig positions, so the offsets and string reads have real locality. A Bloom
  probe is uniformly random by construction, so it is a guaranteed miss with no
  reuse -- it *adds* an access pattern strictly worse than the one it replaces.
- **Gate A already said so.** The tiny-dict filter caps at 8 MB on cache
  residency grounds. Over 96.3 M k-mers even 1 bit/key is 16 MB, so that gate
  would have rejected every configuration here. The prefilter works on a 1.5 M
  k-mer probe panel (2 MB, L2-resident) and does not generalize upward.

The lesson generalizes past Bloom filters: on this path, adding *any*
hash-addressed side structure trades a local access for a random one, and the
existing chain is local enough that the trade loses.

## 6. Sub-1% A/B differences between binaries are not trustworthy

Two separately-built binaries differ in code and data layout, and that shows up
as a **systematic** offset -- stable across runs, so extra repetitions do not
average it away and it reads exactly like a real effect.

Measured directly. Two builds differing only in the sshash SPSS decode, a
function the tiny backend never calls, were compared on the k=23 probe panel
(150M read pairs, `--dict auto` -> tiny, so the differing code cannot execute):

| | mean CPU | pairs lost |
|---|---|---|
| build A | 263.62 s | -- |
| build B | 265.54 s (+0.73%) | 4 of 5 |

An earlier round of the same comparison gave +0.77%, 3 of 4. So the floor for
between-binary comparison on this workload is around **1%**, and it is
directional and repeatable.

Within a single binary, repetition is fine -- the spread is ~0.5% -- so this is
a layout artifact, not run-to-run noise.

Consequences for reading anything else here:

- Effects **above ~3%** are safe to attribute (bucket memoization at -11.6%,
  the k=31 spill fold at -3.2%).
- Effects **near 1%** are not attributable to the code under test without
  building each variant several times with perturbed layout. A tiny-dict decode
  change measured -1.2% on gencode and +1.4% on Flex; both are inside this floor
  and neither number survives. It was reverted (`a0640b8`) on the grounds that
  no benefit was demonstrated, not because a regression was proven.
- The larger input helps precision but **not** this bias: it shrinks run-to-run
  spread while leaving the layout offset untouched.

Use `-t 1` byte-identical RAD output to confirm correctness, and reserve
performance claims for effects several times the floor.

## 7. Environment overrides

| variable | effect |
|---|---|
| `PISCEM_DECODE_SLOTS` | fixes the aggregate decoder slot allocation and disables adaptation; intended for oracle sweeps |
| `PISCEM_RAPIDGZIP_THREADS` | legacy per-decoder-capable-input fixed allocation; disables adaptation |
| `PISCEM_TINY_PREFILTER` | `0` forces the filter off, `1` forces it on |
| `PISCEM_TINY_PREFILTER_MAX_BYTES` | overrides Gate A |
| `PISCEM_TINY_PREFILTER_MIN_MISS` | overrides Gate B's threshold |

Decoder fixed requests must be positive integers. Setting both decoder variables
is an error. Requests above the effective budget are clamped with a warning to
leave one mapping slot, and the requested and applied values are retained in
`map_info.json`.
