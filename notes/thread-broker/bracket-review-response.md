# Response to `bracket-review.md`

The reported t32 bimodality was reproduced and accepted as a release blocker.
The final correction passes the requested eight-repetition t32 gate while
preserving the opposite t8 response shape.

## Evaluation

The review correctly identified that `e78757c` treated an inconclusive first
model point as permission to restore the opening. Because the opening was only a
hint, that made noisy startup evidence capable of reinstating the exact
producer-four answer the bracket was intended to retire.

The review's first proposed correction—require a material, statistically
separated loss to reject the model—was necessary but not sufficient on its own.
One of eight new t32 runs still produced a formally separated false regression.
Applying the proposal literally also broke the opposite t8 case: two t8 runs
were inconclusive between opening four and model one, retained one, and slowed to
18.95--19.23 seconds instead of discovering producer five.

This established that the opening/model comparison is a useful trigger but is
not a stationary population capable of choosing the final fallback on both
surfaces.

## Final design

The implemented policy is:

1. Measure the model point over the configured short opening horizon.
2. If it appears regressed, extend only that rejection attempt to the ordinary
   ratification sample count. The common accept path remains short.
3. Persistent empirical slack/source cap evidence at or below the model target
   independently confirms it and suppresses redundant exploration.
4. Otherwise a regressed or inconclusive result may trigger adjacent points.
5. Each adjacent point must beat the higher measured opening/model rate by more
   than 5% with non-overlapping intervals.
6. If no alternative wins, mapping ends, or the total deadline expires, retain
   the model—not the opening.

The opening is therefore a pivot for candidate placement and a high-water
evidence bar. It is never the fallback. This directly satisfies the review's
central rule that inconclusive evidence cannot restore an unevidenced opening.

The cap discriminator is generic runtime evidence rather than a new modality or
budget constant. t32 establishes a useful producer cap at the model target;
t8 retains headroom and therefore still tests producer five.

## Final real-data gate

Release binary SHA-256:
`78cebc9c499e93ee93e96abe520a786cb9efd18238cfb742e89ec41ed6b03317`.

All runs used the 2,000,000-read / 528,410-mapped scATAC fixture.

### t32, eight repetitions

- Final allocation: producer 1 in 8/8 runs; no producer-four outcome.
- Outcome histogram: `model_selected` 7, `budget_exhausted` 1.
- Mapping wall: 4.526--11.859 s; median 7.808 s.
- Median/reference ratio: 0.907 against the reviewer's 8.61 s pin-one median;
  the required <=1.05 gate passes.
- Bracket wall: median 0.800 s; seven normal paths 0.700--1.001 s; maximum
  4.000111 s.
- The maximum is 111 microseconds beyond the nominal wake deadline due to
  scheduler timing. That run retained producer one rather than restoring four.
- Controller CPU: 0.129--0.363 ms; median 0.173 ms.

Evidence:
`/scratch3/rob/tmp/tb-bracket-cap-gate-final-t32-r{1..8}/`.

### t8 guard, two repetitions

- Final allocation/outcome: producer 5 / `alternative_selected` in 2/2 runs.
- Mapping wall: 14.412--14.839 s; median 14.625 s.
- Bracket wall: 1.052--1.153 s.
- Controller CPU: 0.296--0.371 ms; median 0.334 ms.

Evidence: `/scratch3/rob/tmp/tb-bracket-cap-gate-final-t8-r{1,2}/`.

Every adaptive output and both contemporaneous t8 pin-five controls produced
canonical digest
`293401b75cfd9a9003e3c83b0fe69070b68246f7c3b5a3ae69296674f36bb909`.

## Rejected intermediate evidence

The failed designs were retained rather than overwritten:

- strict rejection without asymmetric confirmation:
  `/scratch3/rob/tmp/tb-bracket-inconclusive-fix-t32-r*`;
- confirmation extension while still restoring the opening:
  `/scratch3/rob/tmp/tb-bracket-confirmed-t{8,32}-r*`;
- model fallback before empirical-cap suppression:
  `/scratch3/rob/tmp/tb-bracket-model-fallback-final-t{8,32}-r*`.

These runs document why simply increasing sample count or treating
inconclusive as an immediate terminal accept did not satisfy both t8 and t32.

## Tests and documentation

Final checks passed:

- `cargo fmt --all -- --check`
- `git diff --check`
- `cargo check --no-default-features`
- strict all-target Clippy with and without `rapidgzip-cpu-accounting`
- 19 controller unit tests
- 45 synthetic control-law tests
- 243 runnable piscem library tests; one fixture test ignored
- three local oracle-boundary tests

New tests cover statistical rejection, cap confirmation, model fallback,
asymmetric confirmation, deadline fallback, and the paired t8/t32 response
shapes.

The design now requires noisy discrete-decision gates to report repetition
count, full range, final-allocation histogram, outcome histogram, and median.
Any cell that has shown bimodality requires at least eight post-fix repetitions,
zero out-of-region outcomes, and a per-run startup-wall bound. A single run is
diagnostic only.

The constant-count pushback is closed as accepted by the review; no additional
constant consolidation was attempted. `.gitignore`, local scripts, manifests,
fixture data, and result trees remain outside the implementation commit.
