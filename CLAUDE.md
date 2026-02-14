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

- **Phase 4: Core Mapping Infrastructure** — Full mapping pipeline from hit collection to RAD output. 54 new tests, 123 total.
  - `src/mapping/hits.rs` — `MappingType`, `HitDirection`, `FragmentEnd`, `SimpleHit`, `SketchHitInfo` trait
  - `src/mapping/sketch_hit_simple.rs` — `SketchHitInfoSimple` (default no-constraint hit tracking)
  - `src/mapping/chain_state.rs` — `ChainState` + `SketchHitInfoChained` (optional structural constraints via SmallVec chains)
  - `src/mapping/filters.rs` — `PoisonState` + `CanonicalKmerIter` + `scan_raw_hits()` (3-phase poison scanning)
  - `src/mapping/cache.rs` — `MappingCache<S: SketchHitInfo>` with `nohash_hasher::BuildNoHashHasher<u32>`
  - `src/mapping/engine.rs` — `map_read<K, S>()` kernel + `collect_mappings_from_hits()` + EC-based ambiguous hit filtering
  - `src/mapping/merge_pairs.rs` — `merge_se_mappings()` paired-end merge (sort + two-pointer + fragment length check)
  - `src/io/rad.rs` — `RadWriter` binary buffer + RAD header/record writers for SC and bulk modes
  - `src/io/fastx.rs` — `FastxSource` wrapping `paraseq` for chunked FASTQ reading (single-end + paired-end)
  - `src/io/threads.rs` — `run_mapping_pipeline()` with crossbeam scoped threads + bounded channel
  - `src/mapping/protocols/mod.rs` — `Protocol` trait + `AlignableReads` + `TechSeqs`
  - New deps: `smallvec`, `nohash-hasher`, `crossbeam`, `paraseq`

- **Phase 5: Protocol Implementations + CLI Wiring** — End-to-end mapping from FASTQ to RAD. 17 new tests, 140 total.
  - `src/io/rad.rs` — REWRITTEN: Correct RAD format (no magic bytes, tag descriptions, file/read/alignment-level tags), `pack_bases_2bit()` (A=0,C=1,G=2,T=3 MSB-first), SC header with `with_position` support, packed BC/UMI (u32/u64), `MappedSecondOrphan` orientation inversion, bulk record with `leftmost_pos` clamping + `u16 frag_len`
  - `src/io/map_info.rs` — `write_map_info()` JSON stats output
  - `src/mapping/map_fragment.rs` — `map_se_fragment()` / `map_pe_fragment()` helpers
  - `src/mapping/protocols/bulk.rs` — `BulkProtocol` with `Protocol` trait impl
  - `src/mapping/protocols/scrna.rs` — `ChromiumProtocol` (V2/V2_5p/V3/V3_5p/V4_3p), `from_name()`, `recover_barcode()`, `barcode_has_n()`, `count_ns()`
  - `src/cli/map_bulk.rs` — Full bulk mapping CLI: load index, write RAD header, `dispatch_on_k!` → `run_bulk_pipeline()`, chunk backpatching, `map_info.json`
  - `src/cli/map_scrna.rs` — Full scRNA mapping CLI: geometry parsing, BC/UMI extraction + N recovery, `--with-position` read length sampling + backpatching
  - New dev-dep: `tempfile`

- **Phase 6: Hardening + Remaining Features** — Unitig-end cache, scATAC protocol, custom geometry, parity harness. 33 new tests, 173 total.
  - `src/mapping/unitig_end_cache.rs` — `UnitigEndCache` with DashMap, capacity-bounded, orientation-aware cache for unitig boundary k-mers shared across threads
  - `src/mapping/streaming_query.rs` — REWRITTEN: Added `cache_end` tracking, `with_cache()` constructor, automatic cache lookup/insert at unitig boundaries
  - `src/mapping/overlap.rs` — Mate overlap detection: `find_overlap()` with dovetail/regular overlap + seed-based alignment with error tolerance
  - `src/mapping/binning.rs` — `BinPos` genome binning for scATAC with overlap regions
  - `src/mapping/protocols/custom.rs` — Hand-written recursive descent parser for custom read geometries (`1{b[16]u[12]x:}2{r:}` format), `CustomProtocol` implementing `Protocol` trait
  - `src/mapping/protocols/scatac.rs` — REWRITTEN: Full `ScatacProtocol` implementing `Protocol` trait
  - `src/io/rad.rs` — Added `write_rad_header_atac()`, `AtacMappingCode` enum, `write_atac_record()` with Tn5 shift
  - `src/cli/map_scatac.rs` — REWRITTEN: Full scATAC CLI with overlap detection, Tn5 shift, genome binning
  - `src/cli/map_scrna.rs` — Updated: custom geometry fallback via `parse_custom_geometry()`, `Box<dyn Protocol>` for polymorphic protocol dispatch
  - `src/verify/index_compare.rs` — REWRITTEN: Real index semantic comparison with `compare_ref_metadata()`
  - `src/verify/rad_compare.rs` — REWRITTEN: RAD file comparison with binary header parser, `compare_rad_files()`, `validate_rad_file()`
  - `src/verify/parity.rs` — REWRITTEN: Parity orchestration with index + RAD comparison, JSON report output
  - `tests/parity_smoke.rs` — REWRITTEN: RAD header roundtrip tests (bulk, SC, ATAC) + ignored integration test
  - New dep: `dashmap = "6"`

- **Phase 7: Poison Table Builder + CanonicalKmer** — `build-poison` CLI command, CanonicalKmer strong type, parity test with poison filtering.
  - `src/mapping/kmer_value.rs` — `CanonicalKmer` newtype wrapping `u64`, with upgrade path docs for k > 31
  - `src/index/build_poison.rs` — Edge-method poison table builder: FASTA scanning with crossbeam threading, `PoisonKmerState` state machine, `scan_sequence_for_poison()` per-sequence scanner
  - `src/cli/poison.rs` — REWRITTEN: Full `build-poison` CLI with dispatch_on_k, index loading, table save
  - `src/index/poison_table.rs` — Updated: `PoisonMap` key and `LabeledPoisonOcc::canonical_kmer` changed from `u64` to `CanonicalKmer`
  - `src/mapping/filters.rs` — Updated: `CanonicalKmerIter` yields `CanonicalKmer`, made `pub(crate)`
  - `src/mapping/unitig_end_cache.rs` — Updated: `DashMap<CanonicalKmer, CachedLookup>`
  - `src/mapping/streaming_query.rs` — Updated: cache key uses `CanonicalKmer`
  - `tests/rad_parity_bulk.rs` — Added: `bulk_pe_rad_parity_with_poison` test (99.68% match rate)
  - Poison table parity: C++ 3,764,601 k-mers / Rust 3,764,549 k-mers (~99.999%)

- **Phase 8: scATAC Parity** — Fix ATAC mapper to match C++ behavior. 100% record-level parity achieved. 7 new tests, 184 total.
  - `src/io/fastx.rs` — Added `ReadTriplet`, `ReadTripletChunk`, `FastxTripleSource` for triple-file FASTQ (R1 + barcode + R2)
  - `src/io/threads.rs` — Added `run_mapping_pipeline_triple()` for triple-file producer-consumer pipeline
  - `src/mapping/hit_searcher.rs` — Added `get_raw_hits_sketch_everykmer()`: queries every k-mer position independently (no contig-walking)
  - `src/mapping/binning.rs` — REWRITTEN: `BinPos` with cumulative per-reference bin IDs matching C++ `bin_pos`, `get_bin_id(tid, pos)` returns `(bin1, bin2)` with overlap region secondary bin
  - `src/mapping/hits.rs` — Added `bin_id: u64` field to `SimpleHit` (default `u64::MAX`)
  - `src/mapping/engine.rs` — Added `map_read_atac<K, S>()`: bin-based hit accumulation with threshold filtering (`ceil(num_valid * thr)`)
  - `src/mapping/merge_pairs.rs` — Added `merge_se_mappings_binned()`: bin-aware PE merge (compatible bins = same or adjacent), `remove_duplicate_hits()` canonicalization, `simple_hit_cmp_bins` comparator, `max_num_hits` post-filter matching C++ `utils.hpp:2018-2027`
  - `src/mapping/map_fragment.rs` — Added `map_se_fragment_atac()` / `map_pe_fragment_atac()` using bin-based mapping + binned merge
  - `src/cli/map_scatac.rs` — REWRITTEN: triple-file input, `--no-poison` defaults true, `--bin-size`/`--bin-overlap`/`--thr` args, every-kmer mode (ignores `--skipping-strategy`), mate overlap → bin-based SE map, no overlap → bin-based PE map + binned merge
  - `examples/atac_mismatch_diag.rs` — Diagnostic tool for categorizing ATAC RAD mismatches
  - `examples/index_stats.rs` — Index statistics dumper (refs, contigs, entries)

### Parity Status

| Mode | Dataset | Mapping Rate | Record-Level Parity |
|------|---------|-------------|-------------------|
| Bulk PE | gencode_pc_v44 (with poison) | 100% match | 99.68% |
| scRNA | SRR12623882 (Chromium V3) | 100% match | 100% |
| scATAC | 5M ATAC reads (hg38 k25) | 100% match (98.33%) | 100% (4,916,721/4,916,721) |

### Next Up

- Performance benchmarking and optimization
- Investigate remaining bulk PE parity gap (0.32% — likely tie-breaking in STRICT/PERMISSIVE mode)

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
      mod.rs                    # build, build_poison, contig_table, eq_classes, poison_table, formats, reference_index, refinfo
      build.rs                  # DONE — End-to-end index build pipeline from cuttlefish output
      build_poison.rs           # DONE — Edge-method poison table builder from decoy FASTA
      contig_table.rs           # DONE — EF offsets + BitFieldVec entries
      refinfo.rs                # DONE — reference names and lengths
      reference_index.rs        # DONE — assembles Dictionary + ContigTable + RefInfo
      eq_classes.rs             # DONE — EC map: tile → EC → (transcript_id, orientation)
      poison_table.rs           # DONE — Poison k-mer table with AHashMap<CanonicalKmer, u64>, build_from_occs, queries
      formats.rs                # ArtifactFormat enum
    cli/
      build.rs                  # DONE — Index build CLI
      poison.rs                 # DONE — build-poison CLI (decoy scanning + save)
      map_bulk.rs               # DONE — Bulk mapping CLI (SE + PE)
      map_scrna.rs              # DONE — scRNA mapping CLI (Chromium protocols, --with-position)
      map_scatac.rs             # DONE — scATAC mapping CLI (overlap detection, Tn5 shift, binning)
    mapping/
      kmer_value.rs             # DONE — CanonicalKmer newtype (u64-backed, upgrade path for k>31)
      unitig_end_cache.rs       # DONE — UnitigEndCache with DashMap<CanonicalKmer>, orientation-aware
      hit_searcher.rs           # DONE — HitSearcher, ReadKmerIter, PERMISSIVE/STRICT modes
      projected_hits.rs         # DONE — RefPos, ProjectedHits<'a>, decode_hit()
      streaming_query.rs        # DONE — PiscemStreamingQuery<'a, K> wrapper + unitig-end cache integration
      hits.rs                   # DONE — MappingType, HitDirection, SimpleHit, SketchHitInfo trait
      sketch_hit_simple.rs      # DONE — SketchHitInfoSimple (no-constraint default)
      chain_state.rs            # DONE — SketchHitInfoChained (optional structural constraints)
      filters.rs                # DONE — PoisonState, scan_raw_hits, CanonicalKmerIter
      cache.rs                  # DONE — MappingCache<S> generic mapping state
      engine.rs                 # DONE — map_read<K,S>() kernel + map_read_atac<K,S>() bin-based kernel
      merge_pairs.rs            # DONE — merge_se_mappings() + merge_se_mappings_binned() bin-aware PE merge
      map_fragment.rs           # DONE — SE/PE helpers + ATAC bin-based variants
      overlap.rs                # DONE — Mate overlap detection (dovetail/regular + seed-based alignment)
      binning.rs                # DONE — BinPos cumulative per-ref binning matching C++ bin_pos
      protocols/
        mod.rs                  # DONE — Protocol trait + AlignableReads + TechSeqs
        bulk.rs                 # DONE — BulkProtocol
        scrna.rs                # DONE — ChromiumProtocol (V2/V2_5p/V3/V3_5p/V4_3p) + barcode recovery
        scatac.rs               # DONE — ScatacProtocol for scATAC-seq
        custom.rs               # DONE — CustomProtocol + geometry parser (recursive descent)
    io/
      rad.rs                    # DONE — RadWriter + RAD headers/records (SC + bulk, with_position)
      fastx.rs                  # DONE — FastxSource + FastxTripleSource (triple-file FASTQ)
      threads.rs                # DONE — run_mapping_pipeline() + run_mapping_pipeline_triple()
      map_info.rs               # DONE — map_info.json writer
    verify/
      mod.rs                    # DONE — Module declarations
      index_compare.rs          # DONE — Index semantic comparison (ref metadata)
      rad_compare.rs            # DONE — RAD file comparison (binary header parser)
      parity.rs                 # DONE — Parity orchestration + JSON report
```

## Terminology Bridge

| sshash-rs term | piscem-cpp term | Meaning |
|---|---|---|
| `string_id` | `contig_id` | Unitig identifier |
| `kmer_id_in_string` | `kmer_id_in_contig` | K-mer position within unitig |
| `num_strings()` | num contigs | Number of unitigs |

## Running Tests

```bash
cargo test              # All 184 tests should pass (12 ignored integration tests)
cargo check             # Should compile clean with no warnings
RUST_LOG=info cargo run # Run with logging
```

## C++ Reference Code

The C++ piscem codebase should be available in `piscem-cpp/` in the current directory for cross-reference. Key files:
- `include/reference_index.hpp` — C++ ReferenceIndex
- `include/basic_contig_table.hpp` — C++ contig table
- `include/hit_searcher.hpp` / `src/hit_searcher.cpp` — Hit collection (~1400 lines)
- `include/mapping/utils.hpp` — `map_read()` kernel
- `include/projected_hits.hpp` — Projected hits and `decode_hit()`
- `include/streaming_query.hpp` — piscem-level streaming query wrapper
