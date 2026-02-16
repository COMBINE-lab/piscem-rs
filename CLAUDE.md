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
- **Phase 10: Multi-file parallel decompression** — Switched from concatenated single-reader streams (`open_concatenated_readers` + `fastq::Reader`) to paraseq's `Collection` API (`fastx::Collection` with `CollectionType::Paired`/`Single`/`Multi`). Enables parallel decompression across multiple input file sets. Processor trait impls updated from `fastq::RefRecord` to `fastx::RefRecord`.

### Parity Status

| Mode | Dataset | Mapping Rate | Record-Level Parity |
|------|---------|-------------|-------------------|
| Bulk PE | gencode_pc_v44 (no poison) | 100% match (96.46%) | 100% (964,594/964,594) |
| Bulk PE | gencode_pc_v44 (with poison) | 100% match | 100% (961,505/961,505) |
| Bulk PE (strict) | gencode_pc_v44 | 100% match | 100% |
| Bulk SE | gencode_pc_v44 | 100% match | 83.65% (tie-breaking differences expected) |
| scRNA | SRR12623882 (Chromium V3) | 100% match | 100% |
| scRNA | PBMC 1k v3 (33.4M reads) | 100% match (86.64%) | 100% (28,968,858/28,968,858) |
| scATAC | 5M ATAC reads (hg38 k25) | 100% match (98.33%) | 100% (4,916,721/4,916,721) |

### Performance Status

Rust is **faster than C++** across both bulk and scRNA workloads (Apple Silicon M2 Max):

**Bulk PE** (1M reads, gencode v44):

| Threads | C++ | Rust | Ratio |
|--------:|----:|-----:|------:|
| 1 | 14.3s | 14.0s | 0.98x |
| 4 | 3.9s | 3.8s | 0.96x |
| 8 | 3.3s | 2.4s | 0.71x |

**scRNA** (PBMC 1k v3, 33.4M reads, Chromium V3, gencode v44, 237K refs):

| Platform | Threads | C++ | Rust | Ratio |
|----------|--------:|----:|-----:|------:|
| Apple Silicon M2 Max | 8 | 114s | 111s | 0.97x |
| x86-64 Linux | 8 | 55s | 47s | 0.85x |

Mapping counts are identical: 28,968,858 / 33,436,697 (86.64%) for both implementations.

Key optimizations applied:
- **AHashMap for hit_map**: Replaced `nohash-hasher` (identity hash) which caused pathological SwissTable H2 collisions with sequential transcript IDs (~38% regression on scRNA with 237K refs). `AHashMap` properly distributes hash bits for SwissTable SIMD probing.
- **AHashSet for observed_ecs**: Replaced standard `HashSet<u64>` (SipHash) with `AHashSet<u64>` matching C++ `ankerl::unordered_dense::set` performance.
- **rapidhash in sshash-rs**: Replaced ahash for MPHF and minimizer hashing. ahash switches algorithm when AES-NI is available (via `target-cpu=native`), silently breaking serialized indices. rapidhash is CPU-feature independent.
- **Optional UnitigEndCache**: Only scATAC uses the cache; bulk and scRNA pass `None`, avoiding DashMap overhead. This was the primary source of the x86-64 performance gap.
- **LocatedHit**: Eliminated double `locate_with_end` Elias-Fano successor queries in dictionary lookups
- **from_ascii_unchecked**: Eliminated `Kmer::from_str` string round-trips (~15% of worker thread time), changed streaming query API from `&str` to `&[u8]`
- **Paraseq native processing**: Zero-copy read access, per-thread stat accumulation (reduced atomic contention at high thread counts)
- **Paraseq Collection**: Multi-file parallel decompression via `fastx::Collection` — threads distributed across reader groups when multiple file sets provided (e.g., `-1 a.fq.gz,b.fq.gz`). No regression for single-file case.

### Next Up

- Benchmark script: `/tmp/bench_piscem.sh` (reads: `test_data/sim_1M_{1,2}.fq.gz`, Rust index: `test_data/gencode_pc_v44_index_rust/`, C++ index: `test_data/gencode_pc_v44_index_cpp/`)
- SC perf test data: `test_data/perf_test/pbmc_1k_v3_S1_L001_R{1,2}_001.fastq.gz`, indices at `test_data/perf_test/{cpp_index,rust_index}`

## Key Design Decisions

1. **No binary compatibility with C++ index format** — Rust has its own serialization format, only semantic equivalence required
2. **Size efficiency matters** — serialized indices should be similar size to C++
3. **Default mapping strategy**: `get_raw_hits_sketch` with PERMISSIVE mode, no structural constraints initially
4. **C++ global mutable state → Rust struct fields**: `ref_shift`/`pos_mask` stored in `EntryEncoding`, passed by reference (not global)
5. **sshash-rs const-generic K**: Use `dispatch_on_k!(k, K => { ... })` at mapping entry point
6. **rapidhash for index hashing**: sshash-rs uses rapidhash (not ahash) for MPHF and minimizer hashing — CPU-feature independent, indices portable across `target-cpu=native`. ahash is still used in piscem-rs for ephemeral hot-path HashMaps (hit_map, observed_ecs) where portability doesn't matter.
7. **UnitigEndCache is scATAC-only**: Bulk and scRNA pass `None` to avoid DashMap overhead
6. **Succinct data structure crates**: `sux` 0.12 git main (with epserde feature), NOT `cseq` or `sucds`
7. **Test data**: Pre-built C++ indices expected in `test_data/` directory
8. **libradicl**: Use git dependency to `develop` branch for RAD comparison
9. **sshash-rs dependency**: Git dependency (`branch = "main"`) in Cargo.toml. Local checkout at `./sshash-rs/` is gitignored and used automatically by Cargo for development.

## Threading Architecture

The mapping pipeline uses **paraseq's `Collection` API** for parallel I/O across multiple input files:

```
paraseq Collection (Vec<fastx::Reader>, CollectionType)
  └─ distributes threads across reader groups (auto: total_threads / num_groups)
     └─ each group: mutex-guarded I/O → fill RecordSet → process batch
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
    let mut readers = Vec::new();
    for (r1, r2) in read1_paths.iter().zip(read2_paths.iter()) {
        readers.push(paraseq::fastx::Reader::new(open_with_decompression(r1)?)?);
        readers.push(paraseq::fastx::Reader::new(open_with_decompression(r2)?)?);
    }
    let collection = Collection::new(readers, CollectionType::Paired)?;
    collection.process_parallel_paired(&mut processor, num_threads, None)?;
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
use paraseq::parallel::{ParallelProcessor, PairedParallelProcessor, MultiParallelProcessor};
use paraseq::fastx::{Collection, CollectionType};
// Processors impl traits for fastx::RefRecord (not fastq::RefRecord) — Collection uses fastx
// impl<'a, 'r, const K: usize> PairedParallelProcessor<paraseq::fastx::RefRecord<'r>> for MyProcessor<'a, K>
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
      fastx.rs                  # open_with_decompression(), Collection/CollectionType re-exports, MultiReader
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
