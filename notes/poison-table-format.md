# Poison table formats (`.poison`)

Both formats hold the same content: each poison k-mer (canonical, ACGT
encoded, first base most significant, `< 2^62`) with its distinct
`(unitig_id, unitig_pos)` occurrences. All integers are little-endian.
`PoisonTable::load` reads either one into the same in-memory table, so mapping
speed and memory do not depend on the file format — only load time and file
size do.

## `PPOIS01` — written by default

```text
"PPOIS01\0"
num_kmers:   u64
num_offsets: u64                         = num_kmers + 1
offsets:     num_offsets × u32           CSR boundaries into occs
num_occs:    u64
occs:        num_occs × {unitig_id: u32, unitig_pos: u32}
entries:     num_kmers × {kmer: u64, offset_index: u64}   sorted by kmer
```

Every piscem-rs release reads it. Loading rebuilds the hash table from
`entries`.

## `PPOIS02` — read, not yet written by default

```text
"PPOIS02\0"
num_kmers: u64   num_occs: u64   num_slots: u64   num_overflow: u64
slots:     num_slots × {key: u64, a: u32, b: u32}
overflow:  num_overflow × {unitig_id: u32, unitig_pos: u32}
```

The in-memory table, verbatim:

- `num_slots = capacity_for(num_kmers) = max(n + n/2, n + 1, 8)` (load factor
  ≤ 2/3). A file with any other slot count is rejected.
- `key = u64::MAX` marks an empty slot. Otherwise the low 62 bits are the
  k-mer, and bit 63 is set when it has more than one occurrence.
- Single occurrence: `a = unitig_id`, `b = unitig_pos`. Multiple:
  `overflow[a .. a + b]`, grouped per k-mer.
- **The slot placement is part of the format.** A k-mer's home slot is
  `fastrange((km ^ (km >> 31)) * 0x9E3779B97F4A7C15, num_slots)`, with linear
  probing. Keys are inserted in increasing k-mer order, so the layout — and
  the file — is deterministic. Changing the hash, the load factor or the
  insertion order needs a new magic (or a rebuild at load, as for
  `PPOIS01`).
- On load, every key must be reachable from its home slot and every overflow
  range in bounds, or the file is rejected as corrupt.

Measured on the GRCh38 decoy of a GENCODE v44 protein-coding index (3.30M
k-mers, 3.46M occurrences): 81.5 MB vs. 93.7 MB for `PPOIS01`; poison load
time ~0.05 s vs. ~0.11 s.

## Switching the default

Older releases fail on a `PPOIS02` file ("invalid poison table magic"), so
only switch `PoisonTable::save` to `PPOIS02` in a release that already
requires indexes to be rebuilt (for example the `PCTAB04`/`PECTB04` bump), and
say so in the release notes. Swap `save`/`save_v2` in `index/poison_table.rs`;
the reader needs no change. Keep reading `PPOIS01` for at least one release
afterwards.
