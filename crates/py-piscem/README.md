# piscem

Python bindings for [piscem](https://github.com/COMBINE-lab/piscem-rs), a fast and accurate tool for k-mer-based read mapping against a reference transcriptome or genome. `piscem` lets you load a pre-built index, map individual reads or read pairs, and perform low-level k-mer queries — all from Python, without running full CLI pipelines.

## Installation

```bash
pip install piscem
```

To build from source (requires Rust and [maturin](https://github.com/PyO3/maturin)):

```bash
cd crates/py-piscem
maturin develop --release
```

## Quick start

```python
import piscem

# Load a pre-built piscem index
index = piscem.ReferenceIndex.load("path/to/index_prefix")

print(f"k={index.k}, {index.num_refs} references, {index.num_contigs} unitigs")
```

## Mapping reads

Create a `MappingEngine` from the index, then map reads one at a time. Each call returns a `MappingResult` containing the mapping type and a list of hits.

### Single-end mapping

```python
engine = index.mapping_engine()

result = engine.map_read(b"ACGTACGT...")
if result.is_mapped:
    for hit in result.hits:
        print(f"{hit.ref_name}  pos={hit.pos}  fw={hit.is_fw}  score={hit.score}")
```

### Paired-end mapping

```python
result = engine.map_read_pair(read1_seq, read2_seq)

if result.mapping_type == "mapped_pair":
    hit = result.hits[0]
    print(f"{hit.ref_name}  frag_len={hit.fragment_length}")
```

The `mapping_type` field is one of `"unmapped"`, `"single_mapped"`, `"mapped_pair"`, `"mapped_first_orphan"`, or `"mapped_second_orphan"`.

### Mapping strategy

Two strategies are available, matching the piscem CLI:

```python
engine = index.mapping_engine(strategy="permissive")  # default — faster, skips along unitigs
engine = index.mapping_engine(strategy="strict")       # queries every k-mer position
```

You can also tune the maximum k-mer occurrence threshold and maximum mappings per read:

```python
engine = index.mapping_engine(max_hit_occ=512, max_read_occ=5000)
```

## Virtual colors mode (binned mapping)

For scATAC-seq and other binned genomic mapping workflows, create a virtual-color engine that accumulates hits per genomic bin:

```python
vc_engine = index.vcolor_engine(bin_size=2000, overlap=400, thr=0.7)

result = vc_engine.map_read_pair(r1, r2)
for hit in result.hits:
    print(f"bin={hit.bin_id}  tid={hit.tid}  pos={hit.pos}")
```

## Streaming k-mer queries

For low-level access, a `StreamingQuery` slides a k-mer window across a sequence and resolves each k-mer against the index, returning unitig and reference coordinates. Consecutive k-mers on the same unitig are resolved by extension rather than a full dictionary lookup.

```python
sq = index.streaming_query()

hits = sq.query_sequence(b"ACGTACGT...")
for i, hit in enumerate(hits):
    if hit is not None:
        for rp in hit.ref_positions:
            print(f"kmer {i}: ref {rp.tid} pos {rp.pos} fw={rp.is_fw}")

print(f"{sq.num_searches} full lookups, {sq.num_extensions} extensions")
```

## Index metadata

```python
index.k               # k-mer size
index.m               # minimizer length
index.num_refs         # number of reference sequences
index.num_contigs      # number of unitigs
index.has_ec_table     # equivalence class table loaded?
index.has_poison_table # poison k-mer table loaded?

index.ref_name(0)      # name of the first reference
index.ref_len(0)       # length of the first reference
index.ref_names()      # list of all reference names
index.ref_lengths()    # list of all reference lengths
```

## Building an index

You can build a new index from [cuttlefish](https://github.com/COMBINE-lab/cuttlefish) output directly from Python:

```python
index = piscem.ReferenceIndex.build(
    "path/to/cuttlefish_prefix",   # .cf_seg, .cf_seq, .json
    "path/to/output_prefix",       # output index files
    k=31,
    m=19,
    threads=8,
)
```

To build an index with a poison k-mer table (for filtering spurious mappings near decoy boundaries), pass one or more decoy FASTA files:

```python
index = piscem.ReferenceIndex.build(
    "path/to/cuttlefish_prefix",
    "path/to/output_prefix",
    k=31,
    m=19,
    threads=8,
    decoys=["path/to/decoys.fa.gz"],
)
print(f"Poison table: {index.has_poison_table}")  # True
```

## Thread safety

- **`ReferenceIndex`** is immutable and can be shared freely across threads.
- **`MappingEngine`** and **`StreamingQuery`** hold mutable per-read state — each thread should create its own via `index.mapping_engine()` or `index.streaming_query()`.

```python
from concurrent.futures import ThreadPoolExecutor

def map_batch(reads):
    eng = index.mapping_engine()
    return [eng.map_read(r) for r in reads]

with ThreadPoolExecutor(max_workers=4) as pool:
    results = list(pool.map(map_batch, read_batches))
```

## API reference

| Class | Description |
|-------|-------------|
| `ReferenceIndex` | Load/build indices, access metadata, create engines |
| `MappingEngine` | Map individual reads or read pairs |
| `MappingResult` | Mapping output: type + list of hits |
| `MappingHit` | Single hit: reference ID, position, orientation, score, fragment info |
| `StreamingQuery` | Low-level sliding-window k-mer queries |
| `KmerHit` | K-mer lookup result: unitig coordinates + reference positions |
| `RefPos` | Position on a reference (tid, pos, orientation) |

## License

BSD 3-Clause
