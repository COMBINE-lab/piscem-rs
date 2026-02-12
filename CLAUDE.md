# piscem-rs Development Context

## Project Overview

This is a Rust port of the C++ `piscem` bioinformatics tool for k-mer-based read mapping. The Rust implementation (`piscem-rs`) must produce **semantically equivalent** outputs to the C++ version (not byte-identical). It depends on `sshash-rs` (local path at `./sshash-rs/`) for the compressed k-mer dictionary.

The full implementation plan with C++ → Rust type mappings, architectural notes, and phased roadmap is in `implementation_plan.md`. Read it before starting new phases.

## Current Status

### Completed Phases

- **Phase 0**: Project bootstrap with CLI skeleton and parity harness scaffolding
- **Phase 1A: ContigTable** (`src/index/contig_table.rs`) — Full implementation with Elias-Fano offsets (`cseq::elias_fano::Sequence`), packed entries (`sux::bits::BitFieldVec`), `EntryEncoding`, `ContigSpan`/`ContigSpanIter`, `ContigTableBuilder`, serialization (magic `PCTAB01\0`). 6 tests passing.
- **Phase 1B: RefInfo** (`src/index/refinfo.rs`) — Reference metadata (names + lengths) with serialization (magic `PRFINF01`). 6 tests passing.
- **Phase 1C: ReferenceIndex** (`src/index/reference_index.rs`) — Assembles `Dictionary` + `ContigTable` + `RefInfo` with `load(prefix, load_ec)`, `save(prefix)`, `from_parts()`, and accessors. Compiles clean, all 13 tests pass.

### Next Up

- **Phase 1D: EqClassMap** (`src/index/eq_classes.rs`) — Currently a stub. Port `equivalence_class_map` with `tile_ec_ids`, `label_list_offsets`, `label_entries`.
- **Phase 1E: Index build pipeline** — Build index from FASTA/cuttlefish input. Parses `.cf_seq` and `.json` files, walks unitigs, builds contig table.
- **Phase 2**: Poison table
- **Phase 3**: Core mapping kernel (PiscemStreamingQuery, ProjectedHits, HitSearcher, mapping cache)
- **Phase 4**: Protocol support (scRNA → bulk → scATAC)
- **Phase 5**: Hardening and performance

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
      mod.rs                    # contig_table, eq_classes, poison_table, formats, reference_index, refinfo
      contig_table.rs           # DONE — EF offsets + BitFieldVec entries
      refinfo.rs                # DONE — reference names and lengths
      reference_index.rs        # DONE — assembles Dictionary + ContigTable + RefInfo
      eq_classes.rs             # STUB — Phase 1D
      poison_table.rs           # STUB — Phase 2
      formats.rs                # ArtifactFormat enum
    cli/                        # CLI subcommands (scaffolded)
    mapping/                    # Mapping engine (scaffolded, not yet implemented)
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
cargo test              # All 13 tests should pass
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
