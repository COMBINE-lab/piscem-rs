# Changelog

Notable changes to `piscem-rs`. Versions follow [Semantic Versioning](https://semver.org);
pre-1.0, a minor bump may carry breaking changes.

## [0.9.0] — unreleased

Supersedes 0.7.0, 0.7.1 and 0.8.0, all of which were **yanked**. Those releases
selected the parallel gzip decoder but never actually decoded in parallel: the
verdict was computed and then had no effect. Mapping output was byte-identical
throughout, so it was a performance failure rather than a correctness one — but it
shipped three times because no test asserted that the decision changed behaviour.
Users on any 0.7.x or 0.8.0 should move to 0.9.0; users on 0.6.4 can upgrade
directly.

### Added

- **Adaptive thread budgeting.** `-t` is now a single execution-slot budget shared
  between mapping and gzip decoding. With `--decoder auto` on seekable gzip input,
  a live thread broker measures each stage's cost while the run proceeds and
  solves for the split rather than guessing it. Measured against swept optima on a
  10x Flex dataset: 5.0% off at `-t 32` and 3.4% off at `-t 64`, settling in 1.4 s
  and 2.6 s. See `docs/threading-and-decompression.md`.
- **Parallel gzip decoding** via a shared decoder pool with one aggregate budget
  across all inputs, replacing per-file thread guesses. On a 150 M read-pair
  dataset this took mapping from 11.28 min to 1.86 min with bit-identical output.
- **`--thread-policy <FILE>`** on all three mapping commands: a JSON document
  overriding thread and decoder policy. Every field is optional; an unknown field
  is an error rather than a silent no-op.
- **`--decoder` for scATAC**, which previously had no way to select a decoder.
- **Machine-readable execution telemetry** in `map_info.json`: the requested and
  effective thread budget, the mapping/decode allocation and how it was chosen,
  the broker's final split and convergence, and consumer cost components.
- **Automatic index-build memory budget.** `build` detects available memory,
  walking the cgroup chain so a container limit is respected rather than the
  host's total.

### Changed

- **`--decoder parallel=N` now honours `-t`.** It previously overwrote the decode
  allocation without reducing the mapping pool, so a pinned request spent both
  sides in full: at `-t 32`, `--decoder parallel=16` ran an average of 41.3
  concurrent threads against a budget of 32. Runs that pin the decoder will now
  use fewer threads than before — correctly — and may take longer as a result.
  `N` remains *per decoder-capable input*; `PISCEM_DECODE_SLOTS=N` sets an
  aggregate and is the better control for reproducible sweeps.
- **The parallel decoder is engaged more conservatively.** It is used only when
  the budget is at least 8 threads per gzip input. Below that the serial decoder
  is faster, because it decodes inline on the mapping threads and so is never
  idle, while a dedicated decode slot cannot map. Measured across 16 cells, the
  serial path won every one where the parallel decoder could not add concurrency,
  by 8–28%. Override with `--thread-policy`.
- **`-t` is capped by available parallelism**, so a cgroup or cpuset limit is
  respected rather than the requested number being taken at face value.
- **Log output defaults to `warn`**, deferring to `RUST_LOG` when it is set.
- The pre-run decoder probe was removed. It predicted the split from a sample and
  landed 7–60% off depending on workload; the broker measures the real run
  instead.

### Fixed

- **RAD output is reproducible.** A run-to-run hash seed made the byte layout of
  equivalent output vary between runs.
- **Non-regular inputs are never handed to the parallel decoder.** A FIFO or
  process substitution (`-r <(zcat …)`) previously risked consuming bytes that
  could not be re-read, and could hang. The downgrade is per file, so a FIFO among
  the inputs no longer costs the regular files their parallel decoder.
- **A failure in the thread broker can no longer fail a run or corrupt its
  output.** A refused resize, an acknowledgement timeout, or a sampler panic
  previously propagated out *after* mapping had completed, unwinding past the RAD
  chunk-count backpatch and leaving a file that parsed cleanly but declared no
  chunks. Broker failures are now advisory, recorded in `map_info.json`.
- Default (non-`rapidgzip`) builds were broken and are restored; both feature
  configurations are covered by CI-equivalent checks.

### Interoperability

- scATAC RAD files name the barcode tag `b`, where C++ piscem used `barcode`;
  reference lengths are written as `u32`, where C++ piscem used `u64`.
  `libradicl` ≥ 0.16 accepts both spellings and the `u32` width, so readers built
  against it — alevin-fry, simpleaf — handle output from either implementation.
  Older readers may not.

### Dependencies

- `rapidgzip-core` 0.3.0, now a released crate rather than a git revision.
- `libradicl` 0.16 → 0.17, `sysinfo` 0.36 → 0.39.
- The `epserde`/`sux`/`value-traits` serialization stack is deliberately held at
  the versions `sshash-lib` requires; they must move together or two copies end up
  either side of the sshash boundary.
- Removed an unused `needletail` entry from the workspace version catalog. The
  version that resolves comes from `sshash-lib` and piscem has no say in it.

## [0.6.4] and earlier

Not tracked in this file. See the commit history.
