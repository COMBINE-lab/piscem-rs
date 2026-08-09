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

**No change.** hashbrown is not uniformly faster — it wins above ~3.6 M keys and
loses below ~1 M — but nothing here justifies a switch:

1. The prefilter already owns the miss path. `lookup_core_bits` records that
   **~79 % of queried k-mers are absent** on probe-panel references and are
   rejected by the blocked Bloom filter in ~12 instructions, before any table is
   consulted. The table choice governs the remaining ~21 %.
2. Where PHast wins (≲ 1 M keys) lookups are L3-resident and already ~4 ns, so
   the win is small in absolute terms.
3. Where PHast loses (≳ 3.6 M keys) it loses by 1.3–1.5×.
4. The memory case is ~16 %, against a 30× build-time cost.

If this is ever revisited, the variant worth measuring is **PHast verifying
against `TinySpss`** rather than a key array: that drops to ~8.3 B/key (a real
2–4× saving) at the cost of one more dependent random access — which is exactly
what already costs PHast the large cases.
