# piscem-rs

A Rust implementation of [piscem](https://github.com/COMBINE-lab/piscem-cpp), a tool for k-mer-based read mapping against a compacted de Bruijn graph index built on the [SSHash](https://github.com/jermp/sshash) data structure.

piscem-rs produces **semantically equivalent** output to the C++ piscem, with 100% record-level parity across all tested modes and datasets.

## Features

- **Index construction** from [cuttlefish](https://github.com/COMBINE-lab/cuttlefish) output (`.cf_seg`, `.cf_seq`, `.json`)
- **Bulk RNA-seq mapping** (single-end and paired-end) to [RAD](https://github.com/COMBINE-lab/libradicl) format
- **Single-cell RNA-seq mapping** with barcode-aware protocols (10x Chromium V2/V3/V4, custom geometries)
- **Single-cell ATAC-seq mapping** with Tn5 shift correction, mate overlap detection, and genome binning
- **Poison k-mer filtering** via decoy-aware index construction
- **Permissive and strict** k-mer skipping strategies with contig-walking acceleration

## Quick start

### Building an index

```bash
# From cuttlefish output
piscem-rs build -i cuttlefish_prefix -o index_prefix -k 31 -m 19

# With equivalence class table
piscem-rs build -i cuttlefish_prefix -o index_prefix -k 31 -m 19 --build-ec-table
```

### Mapping reads

```bash
# Bulk paired-end
piscem-rs map-bulk -i index_prefix -1 reads_1.fq.gz -2 reads_2.fq.gz -o output_dir -t 16

# Single-cell RNA (10x Chromium V3)
piscem-rs map-scrna -i index_prefix -1 reads_1.fq.gz -2 reads_2.fq.gz -o output_dir -t 16 -g chromium_v3

# Single-cell ATAC
piscem-rs map-scatac -i index_prefix -1 R1.fq.gz -b barcode.fq.gz -2 R2.fq.gz -o output_dir -t 16
```

## Building from source

Requires Rust 1.85+.

```bash
git clone --recursive https://github.com/COMBINE-lab/piscem-rs.git
cd piscem-rs
cargo build --release
```

The binary will be at `target/release/piscem-rs`.

## Parity with C++ piscem

piscem-rs is validated against C++ piscem using record-level RAD output comparison:

| Mode | Dataset | Mapping Rate | Record Parity |
|------|---------|:------------:|:-------------:|
| Bulk PE | gencode v44 (1M reads) | 96.46% | 100% |
| Bulk PE + poison | gencode v44 (1M reads) | 96.15% | 100% |
| Bulk PE strict | gencode v44 (1M reads) | 96.46% | 100% |
| scRNA | SRR12623882 (Chromium V3) | — | 100% |
| scATAC | 5M ATAC reads (hg38 k25) | 98.33% | 100% |

## Performance

Single-threaded mapping performance on 1M paired-end reads (gencode v44, Apple Silicon):

| Threads | C++ | Rust | Ratio |
|--------:|----:|-----:|------:|
| 1 | 13.9s | 14.7s | 1.06x |
| 4 | 3.8s | 3.7s | 0.98x |
| 8 | 3.3s | 3.1s | 0.94x |

## Architecture

piscem-rs uses a modular architecture:

- **`sshash-rs`** — Rust port of the [SSHash](https://github.com/jermp/sshash) compressed k-mer dictionary, with streaming query support and [PHast](https://github.com/jermp/pthash) minimal perfect hash functions
- **Index layer** — ContigTable (Elias-Fano + packed entries), RefInfo, EqClassMap, PoisonTable
- **Mapping engine** — Sketch-based hit collection with permissive/strict k-mer skipping, paired-end merge, poison filtering
- **Protocol layer** — Pluggable protocol trait for bulk, scRNA, and scATAC workflows
- **I/O layer** — RAD binary output, chunked FASTQ reading via [paraseq](https://crates.io/crates/paraseq), crossbeam thread pool

## License

MIT
