# piscem-rs Development Context

## Project Overview

This is a Rust port of the C++ `piscem` bioinformatics tool for k-mer-based read mapping. The Rust implementation (`piscem-rs`) must produce **semantically equivalent** outputs to the C++ version (not byte-identical). It depends on `sshash-rs` (local path at `./sshash-rs/`) for the compressed k-mer dictionary.

The full implementation plan with C++ → Rust type mappings, architectural notes, and phased roadmap is in `implementation_plan.md`. Read it before starting new phases.

## Current Status

### Completed Phases

- **Phase 0**: Project bootstrap with CLI skeleton and parity harness scaffolding
- **Phase 1A: ContigTable** (`src/index/contig_table.rs`) — Full implementation with Elias-Fano offsets (`cseq::elias_fano::Sequence`), packed entries (`sux::bits::BitFieldVec`), `EntryEncoding`, `ContigSpan`/`ContigSpanIter`, `ContigTableBuilder`, serialization (magic `PCTAB01\0`). 6 tests passing.
- **Phase 1B: RefInfo** (`src/index/refinfo.rs`) — Reference metadata (names + lengths) with serialization (magic `PRFINF01`). 6 tests passing.
- **Phase 1C: ReferenceIndex** (`src/index/reference_index.rs`) — Assembles `Dictionary` + `ContigTable` + `RefInfo` + optional `EqClassMap` + optional `PoisonTable` with `load(prefix, load_ec, load_poison)`, `save(prefix)`, `from_parts()`, and accessors.
- **Phase 1D: EqClassMap** (`src/index/eq_classes.rs`) — Full implementation with `Orientation` enum, `EcSpan`/`EcSpanIter`, `EqClassMap` (tile → EC → label entries), `EqClassMapBuilder`, serialization (magic `PECTB01\0`). Integrated into `ReferenceIndex` as `Option<EqClassMap>`. 8 tests passing.
- **Phase 1E: Index build pipeline** (`src/index/build.rs`) — Full end-to-end build from cuttlefish output (.cf_seg, .cf_seq, .json). `BuildConfig` + `build_index()` entry point. Parses segments via `sshash_lib::parse_cf_seg()`, builds SSHash dictionary via `DictionaryBuilder`, two-pass .cf_seq parsing (ref info + contig table population), EC table construction from contig entries, short reference handling from JSON. CLI wired up in `src/cli/build.rs`. 10 new tests, 31 total.

- **Phase 2: PoisonTable** (`src/index/poison_table.rs`) — Full implementation with `AHashMap<u64, u64>` (fixed-seed `ahash::RandomState` for deterministic fast lookup), `PoisonOcc`/`LabeledPoisonOcc`, `build_from_occs()`, query methods (`key_exists`, `key_occurs_in_unitig`, `key_occurs_in_unitig_between`), serialization (magic `PPOIS01\0`), JSON stats output. Integrated into `ReferenceIndex` as `Option<PoisonTable>`. CLI stub awaits Phase 3 streaming query engine. 12 new tests, 43 total.

- **Phase 3A+3B: ProjectedHits + PiscemStreamingQuery** — `ProjectedHits<'a>` (`src/mapping/projected_hits.rs`) with `RefPos`, 4-case `decode_hit()` orientation logic, accessors, `resulted_from_open_search` flag. `PiscemStreamingQuery<'a, K>` (`src/mapping/streaming_query.rs`) thin wrapper around sshash-rs `StreamingQueryEngine`. `ReferenceIndex::resolve_lookup()` bridges `LookupResult` → `Option<ProjectedHits>`. 8 new tests, 51 total.

- **Phase 3C: HitSearcher** (`src/mapping/hit_searcher.rs`) — Core k-mer hit collection engine. `SkippingStrategy` enum (Strict/Permissive), `KmerMatchType` enum, `ReadKmerIter` (N-skipping, Clone for save/restore), `HitSearcher` struct with `get_raw_hits_sketch<K>()`. PERMISSIVE mode: skip along contigs using SPSS verification (`check_direct_match`), binary search midpoint recovery, index query fallback. STRICT mode: delegates to `walk_safely_until` which queries every position and extends along contigs via SPSS comparison. Added `Dictionary::kmer_at_pos()` to sshash-rs for SPSS k-mer access. Added `ProjectedHits` setter methods (`set_global_pos`, `set_contig_pos`, `set_contig_orientation`). Unitig-end cache deferred to Phase 5. 18 new tests, 69 total.

### Next Up

- **Phase 4**: Protocol support (scRNA → bulk → scATAC)
- **Phase 5**: Hardening and performance (unitig-end cache, etc.)

## Key Design Decisions

1. **No binary compatibility with C++ index format** — Rust has its own serialization format, only semantic equivalence required
2. **Size efficiency matters** — serialized indices should be similar size to C++
3. **Default mapping strategy**: `get_raw_hits_sketch` with PERMISSIVE mode, no structural constraints initially
4. **C++ global mutable state → Rust struct fields**: `ref_shift`/`pos_mask` stored in `EntryEncoding`, passed by reference (not global)
5. **sshash-rs const-generic K**: Use `dispatch_on_k!(k, K => { ... })` at mapping entry point
6. **Succinct data structure crates**: `sux` (with epserde feature), `cseq`, NOT `sucds`
7. **Test data**: Pre-built C++ indices expected in `test_data/` directory
8. **libradicl**: Use git dependency to `develop` branch for RAD comparison

## Critical Import Patterns

These caused compilation errors previously — remember them:

```rust
// BitFieldVec requires these trait imports for index_value/set_value/len:
use value_traits::slices::{SliceByValue, SliceByValueMut};

// epserde serialization:
use epserde::ser::Serialize;
use epserde::deser::Deserialize;

// CseqSequence size_bytes() requires:
use dyn_size_of::GetSize;
```

## Project Structure

```
piscem-rs/
  Cargo.toml                    # Main crate config
  implementation_plan.md        # Detailed roadmap (READ THIS)
  sshash-rs/                    # Local sshash-rs checkout (path dependency)
    crates/sshash-lib/          # SSHash dictionary library
      src/builder/cf_seg.rs     # Cuttlefish .cf_seg parser (we added this)
  src/
    lib.rs                      # Modules: cli, index, io, mapping, verify
    main.rs                     # Entry point
    index/
      mod.rs                    # build, contig_table, eq_classes, poison_table, formats, reference_index, refinfo
      build.rs                  # DONE — End-to-end index build pipeline from cuttlefish output
      contig_table.rs           # DONE — EF offsets + BitFieldVec entries
      refinfo.rs                # DONE — reference names and lengths
      reference_index.rs        # DONE — assembles Dictionary + ContigTable + RefInfo
      eq_classes.rs             # DONE — EC map: tile → EC → (transcript_id, orientation)
      poison_table.rs           # DONE — Poison k-mer table with AHashMap, build_from_occs, queries
      formats.rs                # ArtifactFormat enum
    cli/                        # CLI subcommands (scaffolded)
    mapping/
      hit_searcher.rs           # DONE — HitSearcher, ReadKmerIter, PERMISSIVE/STRICT modes
      projected_hits.rs         # DONE — RefPos, ProjectedHits<'a>, decode_hit()
      streaming_query.rs        # DONE — PiscemStreamingQuery<'a, K> wrapper
    io/                         # I/O utilities (scaffolded)
    verify/                     # Parity verification (scaffolded)
```

## Terminology Bridge

| sshash-rs term | piscem-cpp term | Meaning |
|---|---|---|
| `string_id` | `contig_id` | Unitig identifier |
| `kmer_id_in_string` | `kmer_id_in_contig` | K-mer position within unitig |
| `num_strings()` | num contigs | Number of unitigs |

## Running Tests

```bash
cargo test              # All 69 tests should pass (1 ignored integration test)
cargo check             # Should compile clean with no warnings
RUST_LOG=info cargo run # Run with logging
```

## C++ Reference Code

The C++ piscem codebase should be available at `../piscem-cpp/` for cross-reference. Key files:
- `include/reference_index.hpp` — C++ ReferenceIndex
- `include/basic_contig_table.hpp` — C++ contig table
- `include/hit_searcher.hpp` / `src/hit_searcher.cpp` — Hit collection (~1400 lines)
- `include/mapping/utils.hpp` — `map_read()` kernel
- `include/projected_hits.hpp` — Projected hits and `decode_hit()`
- `include/streaming_query.hpp` — piscem-level streaming query wrapper
