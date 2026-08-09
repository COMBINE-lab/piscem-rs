# Response to the thread-broker policy-knob review

The review has been evaluated and acted on. The implementation is committed as
`e78757c` (`feat: self-bracket thread broker openings`) and has been
fast-forwarded into the local `dev` branch. Together with the preceding reviewed
t32 correction, local `dev` is two commits ahead of `origin/dev`; it has not been
pushed.

## Changes made

- Removed the performance-derived scATAC producer floor. The actual shared-pool
  safety floor remains one at every budget.
- Replaced the four public nonlinear-search configuration knobs with an opt-in
  `OpeningPolicy::Bracket(OpeningBracketConfig)`.
- Made the scATAC opening a uniform four producer slots at t8, t32 and t64. It
  is explicitly a leaveable hint, not a constraint.
- Removed the `budget <= 8` controller-policy branch.
- Made the startup bracket run only after a stable model/opening disagreement.
  It tests the model answer first. If that answer loses, it restores the opening
  and tests bounded adjacent candidates, first away from and then toward the
  rejected model.
- Required an optional adjacent candidate to be statistically separated and
  more than 5% faster before retaining it.
- Bounded the generic default bracket to three candidate points, 200 ms of
  evidence per point, and a four-second total decision budget.
- Kept the evidence horizon separate from the ordinary resize blackout. Reusing
  the evidence horizon as a blackout added roughly 400 ms to the t8 startup path
  without adding evidence.
- Added public bracket telemetry: measured points, samples, wall time, and
  outcome.
- Ensured the bracket is startup-only and cannot be rearmed by a steady-state
  resurvey. Settled cadence and recurring work are unchanged.
- Ensured every generated candidate respects both configured safety floors.
- Added a mechanical oracle check: a fixed winner is not called an oracle unless
  the swept grid contains a measured point on both sides. Boundary winners now
  fail qualification as unbracketed candidates.

## Measured acceptance results

All real runs used the two-million-read / 528,410-mapped scATAC fixture and
release binary SHA-256
`481b2be0163aa15c9feed1e1ac4f8c31d9f0c309447cf70809e3fa78f7380126`.

| Budget | Mapping wall | Decision | Bracket cost | Controller CPU | Result |
| ---: | ---: | --- | --- | ---: | --- |
| t8 | 13.620, 14.877 s; median 14.248 s | 4 -> 1 -> 4 -> 5; final 5 | 2 points, 16 samples, 1.353--1.453 s | 0.404--0.465 ms | Median/oracle 1.0357 versus the 13.7575 s pin-5 median; the <=5% gate passed. |
| t32 | 7.938 s | 4 -> 1; final 1 | 1 point, 2 samples, 0.800 s | 0.205 ms | Faster than the reviewer's 8.61 s pin-1 median and inside the accepted producer 1--3 region. |
| t64 | 4.142, 4.618 s; median 4.380 s | 4 -> 1; final 1 | 1 point, 2 samples, 0.701--0.800 s | Below the noisy same-binary pin-1 median of 4.925 s. |

At t8, median controller CPU was 0.434 ms, equivalent to 0.00305% of one
core over mapping wall. Every checked adaptive and fixed output retained
canonical digest
`293401b75cfd9a9003e3c83b0fe69070b68246f7c3b5a3ae69296674f36bb909`.

Evidence is retained under
`/scratch3/rob/tmp/tb-opening-bracket-final-t{8,32,64}-*`. The rejected 300 ms
horizon trials remain available as t8 r3/r4; retaining both the passing and
failing repetitions prevented selection of a favorable single run.

## Other modalities and steady-state cost

Bulk SE, bulk PE, scRNA and Flex retain `OpeningPolicy::Fixed`; configuration
tests prove they cannot enter the new bracket. The bracket is never rearmed
after startup, so it adds no settled sampling or administrative work. The
previously measured settled-policy costs remain applicable:

| Policy | Direct controller + sampler CPU | Fraction of one core | Interpretation |
| --- | ---: | ---: | --- |
| Responsive, 25 ms | 4.965 ms | 0.02936% | Higher-frequency regime detection. |
| Sparse-responsive, 5 s | 1.408 ms | 0.00835% | Closest policy to the original low-overhead responsive design. |
| Model-only freeze | 0.401 ms | 0.00197% | Near-free, but unsafe for the non-monotone scATAC response surface. |
| Freeze after full calibration | <=2.582 ms | <=0.0189% over 13.7 s | Finds the response-curve region, then eliminates recurring controller work. |

No new cross-modality real-data matrix was run because those modalities retain
the fixed opening policy and the project had already selected a focused
validation scheme for updates. The normal and comprehensive real matrices
remain available as opt-in pre-release or risk-triggered validation.

## Validation completed

- `cargo fmt --all -- --check`
- `git diff --check`
- `cargo check --no-default-features`
- Strict all-target Clippy with and without `rapidgzip-cpu-accounting`
- 16 controller unit tests
- 44 synthetic control-law tests
- 243 runnable piscem library tests; one fixture test ignored
- Three local oracle-boundary tests
- Real t8/t32/t64 split, telemetry, and canonical-output checks

## Constant-count criterion

The review's category-error finding is fully addressed, but its desired count
of no more than two total scATAC numeric choices was not applied literally.
The surviving values have independent semantics:

- opening four: a leaveable performance hint;
- reader batch 1024: drain and throughput geometry;
- progress publication 64: measurement resolution;
- steady probe 5 s: the independently measured sparse-responsive policy.

Only the opening and sparse policy enter the broker configuration, and none of
these values is a safety floor. Folding batching, measurement resolution, and
settled responsiveness together merely to reduce the count would recreate the
same semantic mistake in another form. Their measurement provenance is recorded
in `design.md` and `completion-ledger.md`.

## Repository state

- Local `dev` head: `e78757c`
- Local `dev`: two commits ahead of `origin/dev`
- Not pushed
- The response, design, completion ledger, oracle scripts, manifests, fixture
  data, and benchmark outputs are intentionally untracked.
- The reviewer's pre-existing `.gitignore` modification remains untouched and
  uncommitted.
