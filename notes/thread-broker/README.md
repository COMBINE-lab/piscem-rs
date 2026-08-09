# Thread broker: design record

Working notes for `crates/thread-broker` and the decode/mapping split it
controls. These are **not** user documentation — for that see
`docs/threading-and-decompression.md` and `docs/thread-broker-release-note.md`.

They are tracked because several shipped constants are only defensible by the
measurements recorded here, and because two of them were previously wrong in ways
that took a full re-measurement to find. A constant whose justification lives in
one person's shell history gets re-guessed rather than re-derived.

## Read in this order

**Start here.** [`design.md`](design.md) — what the controller is, why it solves
for the split rather than searching for it, and the failure of the law it
replaced. §11 lists where it is weakest; §12 says what would falsify it.

### The measurement problem, and the audit

| | |
|---|---|
| [`plan-adaptive-scheduling.md`](plan-adaptive-scheduling.md) | the original validation and measurement strategy |
| [`audit.md`](audit.md) | independent audit; defines the TB-01…TB-12 findings and the A–H gates the rest of these notes refer to |
| [`completion-ledger.md`](completion-ledger.md) | live execution record against that audit — retains failures, including ones later withdrawn |

### The review chain

Each review is followed by its response. Read them in pairs.

| review | response | subject |
|---|---|---|
| [`review-2.md`](review-2.md) | [`review-2-response.md`](review-2-response.md) | post-audit implementation. Found a defect where a broker failure corrupted the output of a *successful* run. Also carries the per-epoch redesign of the TB08 ratification gate |
| [`policy-knobs.md`](policy-knobs.md) | [`policy-knobs-response.md`](policy-knobs-response.md) | scATAC's hard-coded constants, and the proposal for a self-bracketing opening |
| [`bracket-review.md`](bracket-review.md) | [`bracket-review-response.md`](bracket-review-response.md) | review of that implementation; found the `t32` bimodality |

### Investigations

| | |
|---|---|
| [`scatac-t32-failure.md`](scatac-t32-failure.md) | the one real performance failure in a 320-run matrix, traced to a policy constant rather than the control law. §4 carries a correction to an earlier wrong diagnosis |
| [`decoder-policy.md`](decoder-policy.md) | when to use the parallel decoder *at all*. Source of `EngagementPolicy`'s default, and of the reversible-decode proposal |

## Recurring lessons

Three mistakes appear more than once here, and are worth knowing before changing
anything:

1. **Do not encode a performance belief in a safety floor.** `min_producer_slots`
   was set to a "measured good" value and thereby made the controller unable to
   express its own model's answer. See `policy-knobs.md`.
2. **A constant fitted to the edge of a sweep is unproven.** Both times a policy
   number was wrong, it was the smallest point in the grid used to justify it.
   `decoder-policy.md` §4 gives the mechanical check.
3. **Divide CPU seconds by wall time and check against the budget before
   believing any comparison.** Three separate invalid measurements in this
   project would have been caught by that one ratio — see `design.md` §9.

## Related

- Experiment tree for reversible decode:
  `../../../EXPERIMENT-reversible-decode/` (temporary, untracked)
- Validation harnesses: `scripts/thread_broker_*.py`
- Local manifests and raw results are deliberately untracked; harness code is not
