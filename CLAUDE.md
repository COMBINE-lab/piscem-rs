# piscem-rs Development Context

## Project Overview

This is a Rust port of the C++ `piscem` bioinformatics tool for k-mer-based read mapping. The Rust implementation (`piscem-rs`) must produce **semantically equivalent** outputs to the C++ version (not byte-identical). It depends on `sshash-rs` (git dependency on `https://github.com/COMBINE-lab/sshash-rs.git`, branch `main`) for the compressed k-mer dictionary. A local checkout at `./sshash-rs/` is used for development and automatically picked up by Cargo.

The full implementation plan with C++ → Rust type mappings, architectural notes, and phased roadmap is in `implementation_plan.md`. Read it before starting new phases.

## Current Status

### Completed Phases

- **Phase 0**: Project bootstrap with CLI skeleton and parity harness scaffolding
- **Phase 1A–1E: Index data structures + build pipeline** — ContigTable (EF offsets + packed entries), RefInfo, ReferenceIndex, EqClassMap, end-to-end build from cuttlefish output
- **Phase 2: PoisonTable** — `AHashMap<CanonicalKmer, u64>` with fixed-seed ahash, serialization (`PPOIS01\0`), query methods
- **Phase 3: Mapping core** — ProjectedHits, PiscemStreamingQuery (sshash-rs wrapper + unitig-end cache), HitSearcher (PERMISSIVE/STRICT modes)
- **Phase 4: Mapping infrastructure** — `map_read<K, S>()` kernel, MappingCache, SketchHitInfo trait, RadWriter, Protocol trait
- **Phase 5: Protocol implementations + CLI** — Bulk/scRNA/scATAC mapping CLIs, ChromiumProtocol, custom geometry parser
- **Phase 6: Hardening** — UnitigEndCache (DashMap), overlap detection, genome binning, parity harness
- **Phase 7: Poison builder + CanonicalKmer** — `build-poison` CLI, CanonicalKmer newtype
- **Phase 8: scATAC parity** — Triple-file input, every-kmer mode, bin-based merge, 100% record parity
- **Phase 9: Idiomatic paraseq refactor** — Replaced custom crossbeam producer-consumer pipeline with paraseq's native `ParallelProcessor`/`PairedParallelProcessor`/`MultiParallelProcessor` traits. Eliminated intermediate `ReadPair`/`ReadTriplet` owned copies; reads processed in-place from paraseq buffers (zero-copy for single-line FASTQ). Per-thread stats flushed once via `on_thread_complete()` instead of per-chunk atomics.

### Parity Status

| Mode | Dataset | Mapping Rate | Record-Level Parity |
|------|---------|-------------|-------------------|
| Bulk PE | gencode_pc_v44 (no poison) | 100% match (96.46%) | 100% (964,594/964,594) |
| Bulk PE | gencode_pc_v44 (with poison) | 100% match | 100% (961,505/961,505) |
| Bulk PE (strict) | gencode_pc_v44 | 100% match | 100% |
| Bulk SE | gencode_pc_v44 | 100% match | 83.65% (tie-breaking differences expected) |
| scRNA | SRR12623882 (Chromium V3) | 100% match | 100% |
| scATAC | 5M ATAC reads (hg38 k25) | 100% match (98.33%) | 100% (4,916,721/4,916,721) |

### Performance Status

Rust is **faster than C++** on bulk PE mapping (1M reads, gencode v44, Apple Silicon M2 Max):

| Threads | C++ | Rust | Ratio |
|--------:|----:|-----:|------:|
| 1 | 14.4s | 14.7s | 1.02x |
| 4 | 3.9s | 3.8s | 0.96x |
| 8 | 3.3s | 2.4s | 0.71x |

Note: at t=1 paraseq uses a single thread for both I/O and mapping (vs the old dedicated producer thread model). At higher thread counts the mutex-based reader model works well with reduced atomic contention.

Key optimizations already applied:
- **LocatedHit**: Eliminated double `locate_with_end` Elias-Fano successor queries in dictionary lookups
- **from_ascii_unchecked**: Eliminated `Kmer::from_str` string round-trips (~15% of worker thread time), changed streaming query API from `&str` to `&[u8]`
- **Paraseq native processing**: Zero-copy read access, per-thread stat accumulation (reduced atomic contention at high thread counts)

### Next Up

- **Profile remaining hot spots**: Use macOS `sample` or `cargo-instruments` to capture a fresh profile and identify the next optimization targets. From the last profile (pre-`from_ascii_unchecked`), remaining hot areas were:
  - Dictionary lookups (~29% of worker time): bucket dispatch, MPHF evaluation, Elias-Fano successor queries
  - `collect_mappings_from_hits` (~15%): hit deduplication and HashMap operations
  - `merge_se_mappings` + PE merge (~24%): sorting and two-pointer merge of SimpleHit vectors
  - These percentages will have shifted now that `from_str` is eliminated — re-profile to get current breakdown
- Benchmark script: `/tmp/bench_piscem.sh` (reads: `test_data/sim_1M_{1,2}.fq.gz`, Rust index: `test_data/gencode_pc_v44_index_rust/`, C++ index: `test_data/gencode_pc_v44_index_cpp/`)

## Key Design Decisions

1. **No binary compatibility with C++ index format** — Rust has its own serialization format, only semantic equivalence required
2. **Size efficiency matters** — serialized indices should be similar size to C++
3. **Default mapping strategy**: `get_raw_hits_sketch` with PERMISSIVE mode, no structural constraints initially
4. **C++ global mutable state → Rust struct fields**: `ref_shift`/`pos_mask` stored in `EntryEncoding`, passed by reference (not global)
5. **sshash-rs const-generic K**: Use `dispatch_on_k!(k, K => { ... })` at mapping entry point
6. **Succinct data structure crates**: `sux` 0.12 git main (with epserde feature), NOT `cseq` or `sucds`
7. **Test data**: Pre-built C++ indices expected in `test_data/` directory
8. **libradicl**: Use git dependency to `develop` branch for RAD comparison
9. **sshash-rs dependency**: Git dependency (`branch = "main"`) in Cargo.toml. Local checkout at `./sshash-rs/` is gitignored and used automatically by Cargo for development.

## Threading Architecture

The mapping pipeline uses **paraseq's native parallel processing traits** — no custom threading infrastructure:

```
paraseq Reader (mutex-guarded I/O)
  └─ spawns N worker threads via std::thread::scope
     └─ each thread: lock reader → fill RecordSet → unlock → process batch
        └─ Processor struct (one clone per thread, lazy-init state)
           ├─ process_record_pair_batch(): map reads, accumulate RAD output
           ├─ on_batch_complete(): backpatch chunk header, flush to shared file
           └─ on_thread_complete(): flush local stats to atomics (once)
```

**Processor pattern** (`src/mapping/processors.rs`):
- Processor structs hold shared `&'a` references (index, cache, output, stats) + `Option<ThreadState>`
- Custom `Clone` impl: copies reference pointers, sets `state: None`
- `ThreadState` lazily initialized on first batch via `get_or_insert_with()`
- `CommonThreadState<'a, K>`: HitSearcher, PiscemStreamingQuery, MappingCache×3, PoisonState, RadWriter, local counters

**Three processor types**:
- `BulkProcessor` — implements `PairedParallelProcessor` (PE) and `ParallelProcessor` (SE)
- `ScrnaProcessor` — implements `PairedParallelProcessor` with BC/UMI extraction
- `ScatacProcessor` — implements `MultiParallelProcessor` for triple-file (R1 + barcode + R2) input

**CLI pattern** (`src/cli/map_*.rs`):
```rust
dispatch_on_k!(k, K => {
    let mut processor = BulkProcessor::<K>::new(index, end_cache, output, stats, strat);
    let reader1 = paraseq::fastq::Reader::new(open_concatenated_readers(&read1_paths)?);
    let reader2 = paraseq::fastq::Reader::new(open_concatenated_readers(&read2_paths)?);
    reader1.process_parallel_paired(reader2, &mut processor, num_threads)?;
});
```

## Critical Import Patterns

These caused compilation errors previously — remember them:

```rust
// BitFieldVec requires these trait imports for index_value/set_value/len:
use value_traits::slices::{SliceByValue, SliceByValueMut};

// epserde serialization:
use epserde::ser::Serialize;
use epserde::deser::Deserialize;

// Elias-Fano (sux-rs):
use sux::dict::elias_fano::{EliasFanoBuilder, EfSeq, EfSeqDict};
use sux::traits::{IndexedSeq, Succ};

// mem_size() for sux-rs types (sux 0.12 / mem_dbg 0.4):
use mem_dbg::{MemSize, SizeFlags};
// ef.mem_size(SizeFlags::default())

// paraseq parallel processing (explicit lifetime on trait impl):
use paraseq::parallel::{ParallelProcessor, PairedParallelProcessor, MultiParallelProcessor, ParallelReader};
// impl<'a, 'r, const K: usize> PairedParallelProcessor<RefRecord<'r>> for MyProcessor<'a, K>
// (NOT RefRecord<'_> — anonymous lifetimes in impl Trait are unstable)
```

## Project Structure

```
piscem-rs/
  Cargo.toml                    # Main crate config (sshash-lib = git dep)
  implementation_plan.md        # Detailed roadmap (READ THIS)
  sshash-rs/                    # Local sshash-rs checkout (gitignored, auto-used by Cargo)
  src/
    lib.rs                      # Modules: cli, index, io, mapping, verify
    main.rs                     # Entry point
    index/
      mod.rs                    # build, build_poison, contig_table, eq_classes, poison_table, formats, reference_index, refinfo
      build.rs                  # End-to-end index build pipeline from cuttlefish output
      build_poison.rs           # Edge-method poison table builder from decoy FASTA
      contig_table.rs           # EF offsets + BitFieldVec entries
      refinfo.rs                # Reference names and lengths
      reference_index.rs        # Assembles Dictionary + ContigTable + RefInfo
      eq_classes.rs             # EC map: tile → EC → (transcript_id, orientation)
      poison_table.rs           # Poison k-mer table with AHashMap<CanonicalKmer, u64>
      formats.rs                # ArtifactFormat enum
    cli/
      build.rs                  # Index build CLI
      poison.rs                 # build-poison CLI (decoy scanning + save)
      map_bulk.rs               # Bulk mapping CLI — creates BulkProcessor, calls paraseq
      map_scrna.rs              # scRNA mapping CLI — creates ScrnaProcessor, calls paraseq
      map_scatac.rs             # scATAC mapping CLI — creates ScatacProcessor, calls paraseq
    mapping/
      processors.rs             # BulkProcessor, ScrnaProcessor, ScatacProcessor (paraseq trait impls)
      kmer_value.rs             # CanonicalKmer newtype (u64-backed, upgrade path for k>31)
      unitig_end_cache.rs       # UnitigEndCache with DashMap<CanonicalKmer>, orientation-aware
      hit_searcher.rs           # HitSearcher, ReadKmerIter, PERMISSIVE/STRICT modes
      projected_hits.rs         # RefPos, ProjectedHits<'a>, decode_hit()
      streaming_query.rs        # PiscemStreamingQuery<'a, K> wrapper + unitig-end cache integration
      hits.rs                   # MappingType, HitDirection, SimpleHit, SketchHitInfo trait
      sketch_hit_simple.rs      # SketchHitInfoSimple (no-constraint default)
      chain_state.rs            # SketchHitInfoChained (optional structural constraints)
      filters.rs                # PoisonState, scan_raw_hits, CanonicalKmerIter
      cache.rs                  # MappingCache<S> generic mapping state
      engine.rs                 # map_read<K,S>() kernel + map_read_atac<K,S>() bin-based kernel
      merge_pairs.rs            # merge_se_mappings() + merge_se_mappings_binned() bin-aware PE merge
      map_fragment.rs           # SE/PE helpers + ATAC bin-based variants
      overlap.rs                # Mate overlap detection (dovetail/regular + seed-based alignment)
      binning.rs                # BinPos cumulative per-ref binning matching C++ bin_pos
      protocols/
        mod.rs                  # Protocol trait + AlignableReads + TechSeqs
        bulk.rs                 # BulkProtocol
        scrna.rs                # ChromiumProtocol (V2/V2_5p/V3/V3_5p/V4_3p) + barcode recovery
        scatac.rs               # ScatacProtocol for scATAC-seq
        custom.rs               # CustomProtocol + geometry parser (recursive descent)
    io/
      rad.rs                    # RadWriter + RAD headers/records (SC + bulk + ATAC)
      fastx.rs                  # open_concatenated_readers(), MultiReader, paraseq re-exports
      threads.rs                # OutputInfo (mutex RAD file), MappingStats (atomic counters), ThreadConfig
      map_info.rs               # map_info.json writer
    verify/
      mod.rs                    # Module declarations
      index_compare.rs          # Index semantic comparison (ref metadata)
      rad_compare.rs            # RAD file comparison (binary header parser)
      parity.rs                 # Parity orchestration + JSON report
```

## Terminology Bridge

| sshash-rs term | piscem-cpp term | Meaning |
|---|---|---|
| `string_id` | `contig_id` | Unitig identifier |
| `kmer_id_in_string` | `kmer_id_in_contig` | K-mer position within unitig |
| `num_strings()` | num contigs | Number of unitigs |

## Running Tests

```bash
cargo test              # All 179 unit tests should pass
cargo check             # Should compile clean with no warnings
RUST_LOG=info cargo run # Run with logging
```

### Parity Tests (require test data + `--features parity-test`)

```bash
cargo test --features parity-test --release --test rad_parity_bulk -- --ignored --nocapture
cargo test --features parity-test --release --test rad_parity_sc -- --ignored --nocapture
```

## C++ Reference Code

The C++ piscem codebase should be available in `piscem-cpp/` in the current directory for cross-reference. Key files:
- `include/reference_index.hpp` — C++ ReferenceIndex
- `include/basic_contig_table.hpp` — C++ contig table
- `include/hit_searcher.hpp` / `src/hit_searcher.cpp` — Hit collection (~1400 lines)
- `include/mapping/utils.hpp` — `map_read()` kernel
- `include/projected_hits.hpp` — Projected hits and `decode_hit()`
- `include/streaming_query.hpp` — piscem-level streaming query wrapper
