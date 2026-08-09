# Hash table choices

Where hashbrown is, why `eq_classes` uses `HashTable`, and a measured answer to
"would PHast be faster than hashbrown in `TinyDictionary`?".

Harness: [`scripts/hb-vs-phast`](../scripts/hb-vs-phast). Machine for every
number below: AMD EPYC 9575F, **32 MB L3**, single-threaded, `opt-level=3 lto=true
codegen-units=1`, min of 7 reps over 2 M queries.

## 1. Where hashbrown actually is

piscem names it once (`index::eq_classes`). Four copies resolve transitively:

| version | pulled by | on a hot path? |
|---|---|---|
| 0.17.1 | **`tiny-dict`**, `indexmap` (noodles, dev-dep), **and us** | yes — the k-mer probe |
| 0.16.1 | `mem_dbg` ← `epserde`/`sshash-lib` | no, size reporting |
| 0.15.5 | `chumsky` ← `seq_geom_parser` | no, parsed once at startup |
| 0.14.5 | `dashmap` 6 | scATAC unitig-end cache only |

**Our own hot maps are not this crate.** `hash::FixedMap` is
`std::collections::HashMap` with a fixed-seed ahash state; `std`'s table is
hashbrown *vendored into the standard library* and moves with the toolchain, not
with Cargo. Swapping them would gain nothing and would re-permute iteration
order — which `hash.rs` exists to pin, because these maps are iterated straight
into RAD output. Don't.

Taking hashbrown 0.17 directly added **no** new copy: it unifies with the one
`tiny-dict` already resolves.

## 2. Why `eq_classes` uses `HashTable`

`EqClassMapBuilder` interns EC labels. Keyed as `HashMap<Vec<(u32, Orientation)>,
u32>` the map has to own a second copy of every distinct label, because `labels`
already owns them all — the densest structure in the build, doubled.

`HashTable<u32>` stores only the id and takes equality from the caller, so the
label lives in one place and the table costs 4 bytes an entry instead of 32 plus
a duplicated payload. `std` cannot express this: `raw_entry` is unstable and
`HashTable` is not in `std` at all. `raw-entry` and `HashTable` are both in
hashbrown 0.17's default feature set.

### The invariant it buys, and what it costs

Storing ids means the table is only meaningful against `labels`: **every id must
index its own label**. Both reads of an id — equality in `find`, rehashing in
`insert_unique` — go through `labels[id]`. Out of bounds panics; *in bounds but
wrong* silently merges two ECs and hands every tile in one of them the other's
transcript set. Only the second is dangerous, because nothing downstream
detects it.

Three things uphold it, and it is worth being precise about which are real:

1. **Compiler-enforced.** `labels` cannot be mutated while a closure reads it:
   the closures borrow it shared while the table is borrowed mutably. The whole
   "mutated during rehash" class is proof, not vigilance.
2. **Structural.** Ids are minted as `labels.len()`, and `labels` is
   append-only — nothing removes, truncates or reorders it, and `build`
   consumes the builder whole. Any future path that removes or reorders labels
   must rebuild this table rather than patch it.
3. **Checked.** `labels.len() as u32` is the one way to get an in-bounds-but-
   wrong id, by wrapping. Unreachable in practice (~137 GB of labels) but
   asserted, because it is the failure no other mechanism catches. The cast
   predates this change; only the consequence is new.

**`add_tile` pushes before inserting, but that ordering is insurance, not a
fix.** Measured on hashbrown 0.17.1 (100 k inserts, 114 684 hasher invocations
from resizes): `insert_unique` passed the *newly inserted* element to the hasher
**zero** times — it is handed that element's hash and never recomputes it. So
reordering the push after the insert works today. It is an implementation
detail rather than a documented guarantee, and the trap is that such a
reordering would compile and pass every test. Pushing first makes the invariant
hold at every instant instead of only at the points hashbrown happens to look.

`test_dedup_survives_table_growth` forces several resizes; the pre-existing
dedup tests are far too small to ever trigger one.

### End-to-end validation of the change

Run 2026-08-09 against a binary built at `7ea5b84`, the commit immediately
before the rewrite. Script: `VALIDATE-hashbrown/validate.sh` (scratch, not
tracked); index rebuilt from regenerated cuttlefish output with the exact flags
the `piscem` wrapper passes to `piscem_rs::cli::build` — **including
`--build-ec-table`, which is not the default and without which none of this
exercises the changed code**.

| stage | result |
|---|---|
| **0. control** — same binary, index built twice | all 7 artifacts byte-identical. Establishes that byte-comparison is meaningful at all |
| **1. index NEW vs OLD** | all 7 byte-identical, including `index.ectab` (242 784 B, the direct product of the rewrite) and `index.tdct` (37 058 573 B) |
| **2. 10 k fixture, `-t 1`** | `map.rad` byte-identical (1 860 453 B), `unmapped_bc_count.bin` byte-identical, 8 995/10 000 mapped — matching the historical 0.6.4 figure |
| **3. full Flex v2, `-t 64`** | 2 302 290 422 read pairs, **1 994 961 121 mapped (86.65 %) on both**, identical to the last digit |

The full-run count also matches `dyn_full/raw.tsv` from July — a different
binary, index build and `-t` semantics — so this is corroborated across eras,
not merely self-consistent.

Discipline checks: CPU/wall was 54.8 (OLD) and 55.1 (NEW) average threads
against a budget of 64, so neither run overspent `-t`. Wall differed by +1.43 %
and CPU by +2.07 %, both under the ~3 % floor below which this project claims
nothing. The broker was genuinely engaged in both (64 controlled slots, final
split 43/21 vs 41/23, steady phase, no errors, no reverts) — this was not a
serial fallback. The move counts differ (1 vs 5) because the controller widens
its deadband with measured uncertainty, which differed run to run; nothing in
this change touches the broker.

Stage 3 cannot be byte-compared: multi-threaded runs interleave records by
thread arrival, which is inherent to the pipeline. Exact content is covered by
stage 2 at `-t 1`.

## 3. hashbrown vs PHast for `TinyDictionary`

`TinyDictionary` indexes k-mers with
`hashbrown::HashMap<u64, u64, rapidhash::fast::SeedableState>`. sshash instead
uses PHast (`ph-temp`) — so: is hashbrown actually the faster lookup?

### An MPHF is not a dictionary

`phast::Function::get` returns a `usize` for *any* input. For a k-mer that is
not in the reference it returns some other key's slot. To answer "present?" a
PHast dictionary must store the key and compare it. The benchmark therefore
gives PHast `Vec<(key, value)>` indexed by the MPHF slot, interleaved so
verification and payload share a cache line — the layout you would really build.

**sshash gets that verification free and tiny-dict cannot.** sshash re-decodes
the k-mer from the SPSS it must store anyway. `tiny-dict`'s
`lookup_core_bits` returns `(string_id, kmer_id, orientation)` straight out of
the packed value and never touches `TinySpss`, so a PHast port would either add
an 8-byte key array (measured below) or add a third dependent random access into
the SPSS. **PHast being right for sshash does not carry over.**

### Measured

Ratio > 1 means PHast is slower. hashbrown load factor is shown because
`with_capacity` rounds buckets to a power of two, so it swings between 44 % and
87.5 % depending purely on where `n` falls.

| keys | hb load | struct sizes | 100 % hit | 50 % | 0 % hit |
|---:|---:|---|---:|---:|---:|
| 917 504 | 87.5 % | 17 / 14 MB — **both in L3** | **0.93** | **0.79** | **0.34** |
| 3 670 016 | 87.5 % | 68 / 57 MB | 1.30 | 1.29 | 1.37 |
| 7 340 032 | 87.5 % | 136 / 114 MB | 1.33 | 1.32 | 1.47 |
| 1 000 000 | 47.7 % | 34 / 16 MB — **only PHast in L3** | **0.47** | 1.11 | 1.56 |
| 8 000 000 | 47.7 % | 272 / 124 MB | 0.95 | 1.46 | 4.47 |
| 16 000 000 | 47.7 % | 544 / 248 MB | 1.04 | 1.52 | 2.46 |

The crossover is **L3 residency, not load factor**. PHast's structures are
smaller, so it wins whenever it fits in 32 MB and hashbrown does not. Once both
spill, hashbrown wins by ~1.3× on hits: its metadata is 1 byte per bucket and
stays cached, and its dependent chain is shorter than PHast's compressed-seed
decode followed by the slot load.

Two artefacts worth not misreading:

- **The 50 % rows are slower than both pure rows** for both structures. That is
  branch misprediction on the hit/miss decision, not table behaviour.
- **Sparse 0 %-hit flatters hashbrown enormously** (2.5 ns vs 18 ns at 4 M). At
  44–48 % load most probes hit an empty control byte and never load an entry.
  At 87.5 % that advantage collapses to 1.37×.

### Memory: the "< 2 bits/key" headline does not survive

The MPHF really is 2.12–2.15 bits/key. But a *dictionary* still pays 8 bytes of
key plus 8 of value, so PHast-as-dictionary is ~16 B/key against hashbrown's
17 B/bucket — a **16 % saving at dense load**, not 2×. The 2.2× at 8 M keys is
not PHast being small, it is hashbrown having just crossed a power of two.

### Build cost

PHast builds **18–53× slower** (7.3 M keys: 112 ms vs 3.6 s).

## 4. Conclusion

**No change — but on one ground only, and not the one it first appeared.**

### The production reference sits where PHast looked best

The Flex v2 probe panel — the actual reason `TinyDictionary` exists — holds
**1 498 349 distinct canonical k-mers** at k=23 (cuttlefish vertex count,
confirmed by `sum_unitig_len − unitigs×(k−1)` = 2 677 043 − 53 577×22). At that
size, on a 32 MB L3:

| | table | vs L3 |
|---|---:|---|
| hashbrown | 2 097 152 buckets × 17 B = **35.7 MB** | just over |
| PHast | 1.5 M × 16 B + 2.1 bits/key = **24.4 MB** | fits |

The on-disk `.tdct` is 35.3 MB, corroborating the hashbrown figure. This is
almost exactly the `n = 1 000 000` sparse row above — 34 MB vs 16 MB, PHast at
**0.47×** on hits. So the real reference is *not* in the ≳3.6 M regime where
hashbrown wins; it is a few megabytes the wrong side of the cache cliff, which
is the single most favourable configuration PHast had.

An earlier draft of this note implied the size range was unfavourable to PHast.
It is not, and that reasoning should not be reused.

### Why it still does not change

1. **The prefilter owns the miss path.** `lookup_core_bits` records that
   **~79 % of queried k-mers are absent** on probe-panel references and are
   rejected by the blocked Bloom filter in ~12 instructions, before the table is
   consulted. A 2× table win applies to the remaining ~21 % of lookups.
2. **The load factor is luck, not design.** 1.5 M rounds to 2 097 152 buckets —
   71.4 % load. A panel 10 % larger stays at 35.7 MB; one 25 % larger doubles to
   71 MB and loses badly. Building a decision on where the current panel happens
   to fall relative to a power of two is exactly the mistake §"Recurring
   lessons" in `thread-broker/README.md` warns about.
3. **Cost.** ~16 % memory at a dense load, against a 30× build-time cost and a
   new failure mode (a stale MPHF silently returning another k-mer's slot).

The item to revisit is therefore not "PHast vs hashbrown" in the abstract but
**PHast verifying against `TinySpss`**, which is the only variant whose memory
win (~8.3 B/key) is large enough to survive point 2.

If this is ever revisited, the variant worth measuring is **PHast verifying
against `TinySpss`** rather than a key array: that drops to ~8.3 B/key (a real
2–4× saving) at the cost of one more dependent random access — which is exactly
what already costs PHast the large cases.
