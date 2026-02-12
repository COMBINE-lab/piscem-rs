# piscem-rs Implementation Plan

## Scope and Ground Rules

This document translates the goals from `initial_plan.md` into a concrete roadmap for implementing `piscem-rs` in Rust, informed by detailed review of the C++ codebase (`piscem-cpp/`) and the Rust `sshash-rs` library.

### Confirmed constraints

1. **Parity target**: Rust outputs must be **semantically equivalent** to `piscem-cpp` outputs (not byte-identical).
2. **Index size target**: Serialized Rust indices should be **very similar in size** to C++ indices (avoid meaningful disk-size bloat).
3. **RAD comparison path**: Use **`libradicl` on the `develop` branch** to read and compare C++ and Rust RAD outputs.
4. **Project structure**: `piscem-rs` should be an **independent top-level project** that depends on `sshash-lib`, not part of the `sshash-rs` workspace.
5. **Dependency source preference**: Prefer local path dependency to the checked-out `sshash-rs`, with optional Git fallback.

---

## High-Level Goals

### G1: Full index parity in Rust
Implement build, serialization, load, and query of the full piscem index:
- SSHash dictionary (`sshash-lib`)
- Inverted tiling index (unitig -> packed occurrence postings)
- Optional equivalence class table
- Optional poison table

### G2: Exact mapping algorithm parity
Implement mapping algorithms exactly as in `piscem-cpp` such that mapping semantics are identical.

### G3: Protocol support priority
Implement in this exact order:
1. Single-cell RNA-seq
2. Bulk RNA-seq
3. Single-cell ATAC-seq

### G4: Performance and storage discipline
Preserve practical performance and keep serialized index sizes near C++.

---

## Proposed Architecture

## Project layout (independent `piscem-rs`)

```text
piscem-rs/
  Cargo.toml
  src/
    lib.rs
    main.rs
    cli/
      mod.rs
      build.rs
      map_scrna.rs
      map_bulk.rs
      map_scatac.rs
      poison.rs
      inspect.rs
    index/
      mod.rs
      reference_index.rs
      contig_table.rs
      eq_classes.rs
      poison_table.rs
      refinfo.rs
      formats.rs
    mapping/
      mod.rs
      engine.rs
      cache.rs
      chain_state.rs
      hits.rs
      merge_pairs.rs
      filters.rs
      streaming_query.rs   # <-- piscem-level streaming query wrapper (NEW)
      hit_searcher.rs       # <-- hit collection algorithms (NEW)
      projected_hits.rs     # <-- projected_hits + decode_hit (NEW)
      protocols/
        mod.rs
        scrna.rs
        bulk.rs
        scatac.rs
    io/
      mod.rs
      fastx.rs
      rad.rs
      threads.rs
    verify/
      mod.rs
      parity.rs
      rad_compare.rs
      index_compare.rs
  tests/
    data/
    parity_index.rs
    parity_query.rs
    parity_mapping_scrna.rs
    parity_mapping_bulk.rs
    parity_mapping_scatac.rs
```

### Core design principles

1. **Algorithm-preserving port first**
   - Prefer direct translation of control flow and state transitions for parity-critical kernels.
   - Optimize after equivalence harness is stable.

2. **Stable semantic contracts at module boundaries**
   - `index::reference_index` exposes read-only query API used by mapping.
   - `mapping::engine` consumes a protocol adapter (scRNA/bulk/scATAC specifics).

3. **Determinism controls for validation**
   - Single-thread deterministic mode for strict regression checks.
   - Multi-thread mode allows output ordering differences only.

---

## C++ → Rust Type Mapping Reference

This section maps key C++ types to their planned Rust equivalents.

### Index types

| C++ type | Rust module | Rust type | Notes |
|---|---|---|---|
| `piscem_dictionary` | `sshash-lib` | `Dictionary` | Already implemented |
| `basic_contig_table` | `index::contig_table` | `ContigTable` | Port EF offsets + compact_vector entries |
| `equivalence_class_map` | `index::eq_classes` | `EqClassMap` | `tile_ec_ids`, `label_list_offsets`, `label_entries` |
| `reference_index` | `index::reference_index` | `ReferenceIndex` | Owns Dictionary + ContigTable + RefInfo + optional EqClassMap |
| `poison_table` | `index::poison_table` | `PoisonTable` | Hash map + offsets + occurrences |
| `ref_sig_info_t` | `index::refinfo` | `RefSigInfo` | Signature metadata |

### Mapping types

| C++ type | Rust module | Rust type | Notes |
|---|---|---|---|
| `projected_hits` | `mapping::projected_hits` | `ProjectedHits` | contigIdx, contigPos, orientation, contigLen, globalPos, k, refRange |
| `simple_hit` | `mapping::hits` | `SimpleHit` | tid, pos, is_fw, num_hits, score |
| `chain_state` | `mapping::chain_state` | `ChainState` | read_start_pos, prev_pos, curr_pos, num_hits, min_distortion |
| `sketch_hit_info` | `mapping::hits` | `SketchHitInfo` | fw/rc chain vectors with structural constraints |
| `sketch_hit_info_no_struct_constraint` | `mapping::hits` | `SketchHitInfoSimple` | Simple counting variant |
| `mapping_cache_info<S,Q>` | `mapping::cache` | `MappingCache` | Per-thread state: hit_map, accepted_hits, query, searcher |
| `hit_searcher` | `mapping::hit_searcher` | `HitSearcher` | 3 variants of k-mer hit collection |
| `piscem::streaming_query<with_cache>` | `mapping::streaming_query` | `PiscemStreamingQuery` | Wraps sshash StreamingQuery + contig table resolution |
| `MappingType` | `mapping::hits` | `MappingType` | Enum: Unmapped, SingleMapped, Orphan1, Orphan2, Pair |
| `poison_state_t` | `mapping::filters` | `PoisonState` | Poison k-mer scan between hit intervals |
| `SkipContext` | `mapping::hit_searcher` | `SkipContext` | Stateful iterator for skip-and-verify hit collection |

### Terminology bridge: sshash-rs → piscem

| sshash-rs term | piscem-cpp term | Meaning |
|---|---|---|
| `string_id` | `contig_id` | Unitig identifier |
| `kmer_id_in_string` | `kmer_id_in_contig` | K-mer position within unitig |
| `string_begin`/`string_end` | contig begin/end | Unitig boundaries in SPSS |
| `string_length()` | contig length | Unitig length |

---

## Key Architectural Notes from Code Review

### 1. The `piscem::streaming_query` layer

**Critical finding**: The C++ `piscem::streaming_query<bool with_cache>` (in `streaming_query.hpp`) is a **piscem-level wrapper** around sshash's streaming query, **not** the same as sshash's own streaming query. It adds:

- **Contig table resolution**: After sshash lookup, resolves contig_id → reference occurrence span via `m_ctg_offsets`/`m_ctg_entries`.
- **Unitig-end k-mer caching**: `boost::concurrent_flat_map<uint64_t, lookup_result>` that caches the lookup result for the last k-mer in each unitig to speed up transitions between unitigs.
- **Direction/extension tracking**: Tracks `m_direction`, `m_remaining_contig_bases`, `m_prev_contig_id` to decide whether to extend (cheap) or do a full lookup.

In Rust, this must be built as `mapping::streaming_query::PiscemStreamingQuery` wrapping `sshash_lib::StreamingQueryEngine`. It takes a reference to the `ContigTable` to resolve spans.

### 2. Hit searcher architecture

The `hit_searcher` (in `hit_searcher.cpp`, ~1400 lines) implements three variants:

- **`get_raw_hits_sketch`** (primary, newer): Uses `SkipContext` for stateful skip-and-verify. Two sub-modes:
  - `STRICT`: Conservative skip with safe-walk fallback.
  - `PERMISSIVE`: Aggressive skip with mid-point verification; on failure, falls back to every-kmer walk for the failing interval.
- **`get_raw_hits_sketch_orig`** (original): Uses `query_kmer_complex()` with safe-skip walk-back.
- **`get_raw_hits_sketch_everykmer`** (exhaustive): Queries every k-mer in the read (no skipping).

The `SkipContext` struct is central — it tracks: read position, contig position, expected k-mer for fast-hit checking, and manages the skip logic. It reads reference k-mers from the SPSS bit vector (`piscem_bv_iterator`). The Rust equivalent will need to use sshash-rs's ability to decode k-mers at specific positions in the SPSS.

### 3. Mapping core (`mapping/utils.hpp`)

The `map_read()` function (~300 lines) is the per-read mapping kernel:
1. Calls `hit_searcher::get_raw_hits_sketch()` to collect raw hits.
2. Optionally runs `poison_state_t::scan_raw_hits()` to check for poison k-mers.
3. Runs `collect_mappings_from_hits` lambda:
   - First pass: builds `hit_map` (tid → SketchHitInfo) from projected hits.
   - If too many hits (> `max_hit_occ`), opens recovery mode with higher threshold.
   - Applies ambiguous-hit EC filtering if EC table is present.
   - Final best-hit selection based on hit count and structural constraints.
4. Populates `accepted_hits` vector of `SimpleHit`.

### 4. Projected hits and `decode_hit()`

`projected_hits::decode_hit(v)` resolves a packed contig table entry `v` into a reference position + orientation. The 4-case orientation logic (contigFW × contigOrientation) is critical for semantic parity.

### 5. Const-generic K propagation

`sshash-rs` uses `dispatch_on_k!(k, K => { ... })` to go from runtime `k` to const generic `K`. In piscem-rs, the `ReferenceIndex` will learn `k` at load time, and all downstream code paths (streaming query, hit searcher, mapping engine) must be invoked inside a `dispatch_on_k!` block. This suggests the main mapping entry point will be generic over `K`, instantiated at load time.

### 6. Global mutable state replacement

C++ uses global-ish patterns:
- `PiscemIndexUtils::ref_shift()` / `PiscemIndexUtils::pos_mask()` — global bit widths for decoding contig entries.
- `CanonicalKmer::k()` — global k-mer size.

In Rust, these should be stored as fields on `ReferenceIndex` (or a shared `IndexParams` struct) and passed through to code that needs them. `ref_shift` and `pos_mask` derive from `ContigTable::ref_len_bits`.

---

## Dependencies

## Required

- `sshash-lib` (local path dependency preferred)
- `needletail` (FASTA/FASTQ)
- `rayon` or `crossbeam` (parallel work scheduling)
- `clap` (CLI)
- `tracing` + `tracing-subscriber` (logging)
- `anyhow` / `thiserror` (error handling)
- `smallvec` (small fixed-capacity vectors in hot paths — needed for `sketch_hit_info` chain vectors)
- `ahash` (fast hash maps where appropriate)
- `memmap2` (efficient large index IO)
- `nohash-hasher` or similar (for integer-keyed maps in hit_map, matching C++ `ankerl::unordered_dense`)
- `sucds` or custom (Elias-Fano and compact_vector implementations for contig table)

## RAD interoperability

- Integrate/bridge to **`libradicl` (`develop`)** for robust RAD read/compare path.
- Implementation options:
  - preferred: Rust-native RAD reader if practical and parity-safe;
  - fallback: FFI bridge to `libradicl` for decode/normalization used by parity tests.

## `sshash-lib` dependency strategy

Primary (local development):

```toml
[dependencies]
sshash-lib = { path = "./sshash-rs/crates/sshash-lib" }
```

Optional CI/repro fallback:

```toml
[dependencies]
sshash-lib = { git = "https://github.com/COMBINE-lab/sshash-rs", package = "sshash-lib" }
```

---

## Detailed Roadmap

## Phase 0 — Project bootstrap + parity harness

**Objective**: Build the infrastructure to verify semantic equivalence early and continuously.

### Deliverables
- Initialize `piscem-rs` crate skeleton and command structure.
- Add local-path `sshash-lib` dependency.
- Build parity harness:
  - invoke C++/Rust runs on same fixtures,
  - normalize outputs,
  - compare semantics (not bytes).
- Add RAD comparison utility path based on `libradicl` develop.

### Exit criteria
- Can run one command that reports parity pass/fail on a toy dataset.

---

## Phase 1 — Index representation and I/O parity

**Objective**: Implement full Rust index loading/saving/query substrate equivalent to C++.

### 1A: ContigTable (core data structure)

Port `basic_contig_table` to Rust:
- `m_ctg_offsets` → Elias-Fano monotone sequence (encodes cumulative posting list boundaries per unitig)
- `m_ctg_entries` → Compact vector of bit-packed entries (each entry = `ref_position | (orientation_bit << ref_len_bits)`)
- `m_ref_len_bits` → `u64` — number of bits for reference position encoding
- `contig_entries(contig_id)` → returns iterator/slice over entries for given unitig

**Decision needed**: Use `sucds` crate for Elias-Fano, or port the C++ `bits::elias_fano` / `bits::compact_vector` directly? (See Q1 in open questions.)

### 1B: RefInfo

Port reference metadata load/save:
- `m_ref_names: Vec<String>` — reference sequence names
- `m_ref_lens: Vec<u64>` — reference sequence lengths
- Serialized as `.refinfo` file (C++ format: binary length-prefixed strings + u64 array)

### 1C: ReferenceIndex

Assemble the full index:
- `Dictionary` (loaded via `sshash-lib`)
- `ContigTable` (loaded from `.ctab` equivalent)
- `RefInfo` (loaded from `.refinfo` equivalent)
- Optional `EqClassMap` (loaded from `.ectab` equivalent)
- `ref_shift` / `pos_mask` computed from `contig_table.ref_len_bits`

Key method: `query(kmer_iter, streaming_query) -> ProjectedHits`
- Takes a k-mer iterator and piscem streaming query
- Returns `ProjectedHits` with contig info + reference range

### 1D: EqClassMap (optional component)

Port `equivalence_class_map`:
- `tile_ec_ids` → compact vector (unitig tile → EC id)
- `label_list_offsets` → Elias-Fano (EC id → offset into label entries)
- `label_entries` → compact vector (reference IDs in each EC)
- Methods: `entries_for_ec(ec_id)`, `entries_for_tile(tile_id)`, `ec_for_tile(tile_id)`

### 1E: Index build from FASTA

Port `build_contig_table.cpp` logic:
- Walk all unitigs from Dictionary, resolve reference occurrences, build packed posting lists
- Construct Elias-Fano offsets and compact vector entries
- Serialize to piscem-rs native format

### Exit criteria
- Rust-built index passes semantic index comparison against C++.
- On-disk index sizes are in expected range (no major inflation).

---

## Phase 2 — Optional poison table parity

**Objective**: Implement poison table generation and query semantics equivalent to C++ optional behavior.

### Implementation targets

Port `poison_table` from `poison_table.hpp`:
- Hash map: canonical k-mer → offset into occurrence array
- Offset array + `poison_occ_t` vector (unitig_id, begin_offset, end_offset entries)
- Query methods: `key_exists()`, `key_occurs_in_unitig()`, `key_occurs_in_unitig_between()`
- Build: `build_from_occs()` — processes raw poison k-mer occurrences
- Serialize: `.poison` file

### Exit criteria
- Poison-enabled runs show semantic parity with C++ reference.

---

## Phase 3 — Core mapping kernel parity

**Objective**: Port exact mapping behavior from C++ core mapping utilities.

### 3A: PiscemStreamingQuery

Port `piscem::streaming_query<with_cache>` as a piscem-rs wrapper around `sshash_lib::StreamingQueryEngine`:
- Holds reference to `ContigTable` for span resolution
- Implements `query_lookup()`: performs sshash lookup, then resolves `contig_span` from contig table
- Tracks direction, remaining contig bases, previous contig ID for extension optimization
- Unitig-end cache: `DashMap<u64, LookupResult>` or similar concurrent map (replaces `boost::concurrent_flat_map`)

### 3B: ProjectedHits and decode_hit

Port `projected_hits` struct and the 4-case `decode_hit(v)` orientation logic:
- Extract orientation bit and reference position from packed entry
- 4 cases: (contigFW × contigOrientation) → (ref_pos, ref_fw)

### 3C: HitSearcher

Port `hit_searcher` (~1400 lines of C++):
- `SkipContext` struct: read_pos, contig_pos, expected_kmer, skip logic
- Primary variant `get_raw_hits_sketch` with STRICT/PERMISSIVE
- Original variant `get_raw_hits_sketch_orig`
- Exhaustive variant `get_raw_hits_sketch_everykmer`
- Left/right raw hit vectors: `Vec<(i32, ProjectedHits)>`

**Port order**: `get_raw_hits_sketch` (PERMISSIVE) first, then STRICT, then orig, then everykmer.

### 3D: Mapping cache and hit accumulation

Port `mapping_cache_info` and `map_read()`:
- `hit_map: HashMap<u32, SketchHitInfo>` (or nohash variant for integer keys)
- `accepted_hits: Vec<SimpleHit>`
- Two-pass collection with occ recovery
- Ambiguous-hit EC filtering
- Best-hit selection

### 3E: Pair merging and filters

Port `merge_se_mappings()`:
- Merge left/right single-end mappings for paired-end reads
- Fragment length computation, concordance checks
- `MappingType` determination (Unmapped/SingleMapped/Orphan1/Orphan2/Pair)

### Exit criteria
- Single-thread mapping outputs semantically identical for controlled fixtures.

---

## Phase 4 — Protocols in priority order

## 4A: scRNA-seq (first)
- Protocol enum: ChromiumV2, V2_5P, V3, V3_5P, V4_3P, Custom
- Barcode and UMI extraction from read sequences based on geometry
- scRNA RAD emission path
- Paired-end mapping dispatch
- Parity test set for multiple protocol variants

## 4B: bulk RNA-seq
- Reuse core engine with bulk-specific options and output semantics.

## 4C: scATAC-seq
- Add scATAC-specific technical-sequence extraction and mapping/reporting behavior.
- `bin_pos` helper for position-based binning.

### Exit criteria
- Per-protocol parity suites passing.

---

## Phase 5 — Performance, storage tuning, and release hardening

**Objective**: Reach production-ready speed/memory behavior without sacrificing equivalence.

### Work items
- Profile hotspots and tune data structures/allocation patterns.
- Validate index size deltas and compactness.
- Add larger regression datasets and CI gates.
- Document reproducibility and benchmark methodology.
- Evaluate unitig-end cache effectiveness and tuning.

### Exit criteria
- Stable parity + acceptable performance + acceptable on-disk size behavior.

---

## Verification Strategy

## Equivalence philosophy

- Compare **meaning**, not bytes.
- For multithreaded outputs, compare as multisets after normalization.

## Index verification

- Compare decoded posting lists for sampled/full unitigs:
  - unitig -> [(orientation, ref_id, pos), ...]
- Compare equivalence-class interpretation.
- Compare reference metadata (`ref names`, `ref lengths`).
- Track size deltas (`rust_size / cpp_size`) as an explicit metric.

## Query verification

- For randomized and fixture k-mer queries, compare:
  - found/not found,
  - projected contig id/offset/orientation,
  - any optional poison behavior.

## Mapping verification

- Compare RAD-derived semantic records through `libradicl` normalization.
- Single-thread: exact semantic identity.
- Multi-thread: order-insensitive equivalence only.

---

## Performance and Size Goals

## Throughput and memory
- Initial target: match C++ within practical range on representative datasets.
- Follow-up target: improve hot-path performance where possible without changing semantics.

## Serialized index footprint
- Keep Rust index disk size close to C++ equivalent.
- Treat persistent, significant size inflation as a release blocker for format/layout tuning.

---

## Implementation TODOs (Actionable)

## Foundation
- [x] Create `piscem-rs` crate scaffold and command skeleton.
- [x] Add local-path `sshash-lib` dependency and compile smoke test.
- [x] Add logging/error conventions and configuration.

## Parity harness
- [ ] Build fixture runner for C++ vs Rust index/query/map.
- [ ] Add normalized comparators for index artifacts.
- [ ] Add RAD semantic comparator using `libradicl` develop path.

## Index
- [ ] Implement ContigTable: Elias-Fano offsets + compact vector entries.
- [ ] Implement RefInfo: ref names/lengths serialization.
- [ ] Implement ReferenceIndex: assemble Dictionary + ContigTable + RefInfo.
- [ ] Implement EqClassMap: tile → EC → label entries.
- [ ] Implement index build pipeline from FASTA input.
- [ ] Implement `decode_hit()` orientation logic (4-case).

## Streaming query layer
- [ ] Implement PiscemStreamingQuery wrapping sshash StreamingQueryEngine.
- [ ] Add contig table span resolution after sshash lookup.
- [ ] Add direction/extension tracking and unitig-end caching.

## Hit searcher
- [ ] Implement SkipContext stateful iterator.
- [ ] Port `get_raw_hits_sketch` (PERMISSIVE mode first).
- [ ] Port `get_raw_hits_sketch` (STRICT mode).
- [ ] Port `get_raw_hits_sketch_orig`.
- [ ] Port `get_raw_hits_sketch_everykmer`.

## Poison
- [ ] Implement poison table builder (edge mode parity first).
- [ ] Integrate poison-aware query/mapping path.
- [ ] Implement `poison_state_t::scan_raw_hits()` equivalent.

## Mapping core
- [ ] Port `map_read()` per-read mapping kernel.
- [ ] Port `sketch_hit_info` with structural constraint chains.
- [ ] Port `sketch_hit_info_no_struct_constraint` simple counting.
- [ ] Port `mapping_cache_info` per-thread cache.
- [ ] Port acceptance filters and tie-breakers exactly.
- [ ] Port `merge_se_mappings()` paired-end merge semantics.

## Protocols
- [ ] Implement scRNA protocol path (barcode + UMI).
- [ ] Implement bulk RNA path.
- [ ] Implement scATAC path.

## Hardening
- [ ] Add comprehensive parity test matrix in CI.
- [ ] Add benchmark suite and disk-size tracking reports.
- [ ] Write user-facing docs for all commands and options.

---

## Open Questions (Awaiting User Input)

1. **Index binary compatibility**: Should Rust load C++-serialized indices (`.sshash`, `.ctab`, `.refinfo`), or build its own format only? Binary compat enables parity testing against the same index but requires matching the exact C++ serialization.

  - No, binary compatibility is not required, only the semantic equivalence of the indices . Though, as noted above, efficiency is paramount so the size of the Rust files should not be much larger than that of the C++ files.

2. **Primary hit searcher variant**: Which `get_raw_hits_sketch` variant/strategy is the default in production use? (Assumed: `get_raw_hits_sketch` with PERMISSIVE.)

  - Yes, this is the default strategy (and the one we want to focus on getting right first)

3. **Test data and C++ binary**: Are toy datasets and a compiled C++ piscem binary available? Or should the parity harness build piscem-cpp from source?

  - There are toy datasets available, though piscem-cpp can be built from source if the executable needs to be run.  In the `test_data` directory, we have 5 imporant things:
    1) a folder `gencode_pc_v44_dbg` that contains the input (segment and sequence) file that is necessary for building the sshash component of the index, and the inverted tiling index, this is the input to the piscem-cpp `build` command
    2) a folder `gencode_pc_v44_index_nopoison`, which are the files generated by running the `build` command without poison/decoy sequences
    3) a folder `gencode_pc_v44_index_with_poison`, which are the files generate dby running the `build` command, and then running `build-poison-table` on the constructed index (and giving `GRCh38.primary_assembly.genome.fa.gz` as the decoy sequence)
    4) the file `GRCh38.primary_assembly.genome.fa.gz` used as the poison/decoy sequence for 2 above
    5) the file `gencode.v49.pc_transcripts.fa.gz`, which we should not need directly, but which is the source file from which the de bruijn graph in 1) was genereated.

4. **`libradicl` integration**: How should we depend on `libradicl` develop? Git dependency? Local path?

  - we should depend on it as a git dependency (to the develop branch). If it is decided we need to modify it for some purpose, I will pull it in locally and we can then depend on and change the local version (whose changes I will push back upstream)

5. **Structural constraints default**: Is `sketch_hit_info` (structural constraints enabled) or `sketch_hit_info_no_struct_constraint` the default? Should we port both from the start?

  - without structural constraints is the default. We will eventually want both, but we should start with no structural constraints

6. **Custom geometry parsing**: How important is `CUSTOM` protocol geometry for the initial implementation?

  - We will eventually want this, but we can delay it until after other features are done

7. **Elias-Fano / compact_vector implementation**: Use `sucds` crate, or port the C++ `bits::` implementations directly? The C++ uses specific template instantiations (`elias_fano<false, false>`, `compact_vector`).

  - Neither; `sshash-rs` already pulls in `sux-rs` and `cseq` and we can feel comfortable relying on the latest version of either (or both) of these. In general, we should first prefer a solution from `sux-rs` (https://github.com/vigna/sux-rs), then anything from `bsuccinct-rs` (https://github.com/beling/bsuccinct-rs).
  
---

## Risk Register and Mitigations

1. **Semantic drift in nuanced mapping logic**
   - Mitigation: golden parity tests at each stage, deterministic single-thread mode.

2. **Index size inflation in Rust serialization**
   - Mitigation: explicit size KPI tracking, bit-packing and succinct structures from start.

3. **RAD parsing mismatch**
   - Mitigation: use `libradicl` develop as canonical decode path for comparison.

4. **Dependency friction between local and CI setups**
   - Mitigation: support both local path and Git dependency configuration for `sshash-lib`.

5. **Const-generic K propagation complexity** (NEW)
   - `sshash-rs` requires compile-time `K` via `dispatch_on_k!`. All downstream code (streaming query, hit searcher, mapping engine) must be generic over `K` or invoked inside a dispatch block.
   - Mitigation: establish the `dispatch_on_k!` boundary early (at `ReferenceIndex::load` or mapping entry point) and keep inner code K-generic.

6. **piscem::streaming_query layer mismatch** (NEW)
   - The sshash-rs `StreamingQuery` does NOT include contig table resolution or unitig-end caching. This is a separate layer that must be built in piscem-rs.
   - Mitigation: clearly separate sshash-level query from piscem-level query in the module structure.

7. **Global mutable state patterns in C++** (NEW)
   - C++ uses `PiscemIndexUtils::ref_shift()`, `CanonicalKmer::k()` as global state. Rust must thread these through as struct fields.
   - Mitigation: define `IndexParams { k, ref_shift, pos_mask }` early and pass by reference.

---

## Immediate Next Steps

1. ~~Create initial `Cargo.toml` + module skeleton for `piscem-rs`.~~ (Done)
2. ~~Wire local-path `sshash-lib` dependency and add build command.~~ (Done)
3. Resolve open questions above (especially Q1: binary compatibility).
4. Implement Phase 1A: `ContigTable` with Elias-Fano offsets and compact vector entries.
5. Implement Phase 1B: `RefInfo` load/save.
6. Implement Phase 1C: `ReferenceIndex` assembling Dictionary + ContigTable + RefInfo.
7. Add first parity test: load same reference, compare decoded posting lists.
