# Poison table formats (`.poison`)

Both formats hold the same content: each poison k-mer (canonical, ACGT
encoded, first base most significant, `< 2^62`) with its distinct
`(unitig_id, unitig_pos)` occurrences. All integers are little-endian.
`PoisonTable::load` reads either one into the same in-memory table, so mapping
speed and memory do not depend on the file format — only load time and file
size do.

## `PPOIS02` — written by default since 0.11.0

```text
"PPOIS02\0"
num_kmers: u64   num_occs: u64   num_slots: u64   num_overflow: u64
seed:      u64
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
  `fastrange(fmix64(km ^ seed), num_slots)` (murmur3's 64-bit finalizer), with
  linear probing, and keys are inserted in increasing k-mer order. Changing
  the hash function, the load factor or the insertion order needs a new magic
  (or a rebuild at load, as for `PPOIS01`). The seed is not fixed; see below.

Measured on the GRCh38 decoy of a GENCODE v44 protein-coding index (3.76M
k-mers, 3.95M occurrences): 93.1 MB vs. 106.9 MB for `PPOIS01`; poison load
time 0.08 s vs. 0.14 s.

### Why the hash cannot make a table slow

The hash is fixed by the format, but a hash that suits one key set can still
cluster another. Linear probing at load 2/3 with a random hash displaces a key
by 1 slot on average (2 probes per hit, 5 per miss), and the longest
displacement at 3.7M keys is 100–160. Real and structured k-mer sets measure
exactly that (genome decoys, consecutive genome k-mers at k = 19/23/31,
microsatellites, homopolymer runs, arithmetic progressions). But a key set can
be built that sends every key to one slot under any one fixed hash, and then
build, load and every query are quadratic. So:

- **Build.** `PoisonTable::layout` tries the seeds `seed_for(0), seed_for(1),
  …` in order and keeps the first whose layout is within a probe budget: total
  displacement ≤ `3n + 1024` and longest displacement ≤ `64 + 16·⌈log2 n⌉`
  (≈ 3× the random-hash mean, and ≈ 2.5× the longest displacement seen over
  6,200 seeds at sizes 40 to 3.7M — none of which exceeded the budget). An
  attempt is abandoned as soon as it exceeds the budget, so every attempt is
  O(n). The sequence is fixed, so a key set always gets the same seed and the
  same file (in practice `seed_for(0)`).
- **Load.** One linear pass checks that every key is reachable from its home
  slot (otherwise the file is rejected as corrupt) and that the layout is
  within the same budget. A valid file that is not — which `save` never
  writes — is laid out again in memory under the next good seed instead of
  being used as is.

So the worst a key set can do is cost an extra O(n) attempt. This replaces
the `(km ^ (km >> 31)) * 0x9E3779B97F4A7C15` hash of the unreleased first
draft of this format, which had no seed and no budget. That draft's file
layout (no `seed` field) was never written by a release.

## `PPOIS01` — still read; written by `save_v1`

```text
"PPOIS01\0"
num_kmers:   u64
num_offsets: u64                         = num_kmers + 1
offsets:     num_offsets × u32           CSR boundaries into occs
num_occs:    u64
occs:        num_occs × {unitig_id: u32, unitig_pos: u32}
entries:     num_kmers × {kmer: u64, offset_index: u64}   sorted by kmer
```

Written by piscem-rs before 0.11.0 (and C++ piscem). Loading rebuilds the
table from `entries`, through the same seeded layout. `PoisonTable::save_v1`
still writes it, byte-identical to earlier releases' files; nothing in the
CLI does.

## Compatibility

Releases before 0.11.0 fail on a `PPOIS02` file ("invalid poison table
magic"). 0.11.0 already requires indexes to be rebuilt (`PCTAB04`/`PECTB04`,
sshash-lib 0.8), so `build`/`build-poison` write `PPOIS02` from 0.11.0 on.
Keep reading `PPOIS01` for at least one release afterwards.
