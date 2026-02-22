//! Python bindings for the piscem k-mer mapping engine.
//!
//! Exposes the following Python classes:
//! - `ReferenceIndex` — load/build indices, metadata, factory methods
//! - `MappingEngine` — per-read SE/PE mapping (standard and virtual-color modes)
//! - `MappingResult` — mapping output with hits
//! - `MappingHit` — single hit on a reference
//! - `StreamingQuery` — low-level k-mer queries
//! - `KmerHit` — result of a k-mer lookup
//! - `RefPos` — position on a reference

use std::path::Path;
use std::sync::Arc;

use pyo3::exceptions::{PyIOError, PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use sshash_lib::{dispatch_on_k, Kmer, KmerBits, LookupResult, StreamingQuery};

use piscem_rs::index::reference_index::ReferenceIndex;
use piscem_rs::mapping::binning::BinPos;
use piscem_rs::mapping::cache::MappingCache;
use piscem_rs::mapping::filters::PoisonState;
use piscem_rs::mapping::hit_searcher::{HitSearcher, SkippingStrategy};
use piscem_rs::mapping::hits::{MappingType, INVALID_FRAG_LEN, INVALID_MATE_POS};
use piscem_rs::mapping::map_fragment::{
    map_pe_fragment, map_pe_fragment_atac, map_se_fragment, map_se_fragment_atac,
};
use piscem_rs::mapping::sketch_hit_simple::SketchHitInfoSimple;
use piscem_rs::mapping::streaming_query::PiscemStreamingQuery;

// ─── MappingOpts ──────────────────────────────────────────────────────────────

/// Shared mapping configuration.
#[derive(Clone)]
struct MappingOpts {
    max_hit_occ: usize,
    max_read_occ: usize,
    max_ec_card: u32,
}

impl Default for MappingOpts {
    fn default() -> Self {
        Self {
            max_hit_occ: 256,
            max_read_occ: 2500,
            max_ec_card: 4096,
        }
    }
}

// ─── DynIndex ─────────────────────────────────────────────────────────────────
//
// Trait-erased dispatch on const-generic K. Dispatched once at load time.

trait DynIndex: Send + Sync {
    fn make_mapping_engine(
        &self,
        strat: SkippingStrategy,
        opts: MappingOpts,
    ) -> Box<dyn DynMappingEngine>;

    fn make_vcolor_engine(
        &self,
        opts: MappingOpts,
        bin_size: u64,
        overlap: u64,
        thr: f32,
    ) -> Box<dyn DynMappingEngine>;

    fn make_streaming_query(&self) -> Box<dyn DynStreamingQuery>;

    fn k(&self) -> usize;
    fn m(&self) -> usize;
    fn num_refs(&self) -> usize;
    fn num_contigs(&self) -> usize;
    fn ref_name(&self, i: usize) -> &str;
    fn ref_len(&self, i: usize) -> u64;
    fn has_ec_table(&self) -> bool;
    fn has_poison_table(&self) -> bool;
    fn save(&self, prefix: &str) -> Result<(), String>;
}

struct ConcreteIndex<const K: usize>
where
    Kmer<K>: KmerBits,
{
    inner: Arc<ReferenceIndex>,
}

impl<const K: usize> DynIndex for ConcreteIndex<K>
where
    Kmer<K>: KmerBits,
{
    fn make_mapping_engine(
        &self,
        strat: SkippingStrategy,
        opts: MappingOpts,
    ) -> Box<dyn DynMappingEngine> {
        Box::new(ConcreteMappingEngine::<K> {
            index: Arc::clone(&self.inner),
            strat,
            opts,
            binning: None,
            cache_left: MappingCache::new(K),
            cache_right: MappingCache::new(K),
            cache_out: MappingCache::new(K),
        })
    }

    fn make_vcolor_engine(
        &self,
        opts: MappingOpts,
        bin_size: u64,
        overlap: u64,
        thr: f32,
    ) -> Box<dyn DynMappingEngine> {
        let binning = BinPos::new(&self.inner, bin_size, overlap, thr);
        Box::new(ConcreteMappingEngine::<K> {
            index: Arc::clone(&self.inner),
            strat: SkippingStrategy::Permissive,
            opts,
            binning: Some(binning),
            cache_left: MappingCache::new(K),
            cache_right: MappingCache::new(K),
            cache_out: MappingCache::new(K),
        })
    }

    fn make_streaming_query(&self) -> Box<dyn DynStreamingQuery> {
        let k = self.inner.k();
        let m = self.inner.m();
        let canonical = self.inner.dict().canonical();
        Box::new(ConcreteStreamingQuery::<K> {
            query: StreamingQuery::new(k, m, canonical),
            index: Arc::clone(&self.inner),
            prev_query_pos: i32::MIN,
        })
    }

    fn k(&self) -> usize {
        self.inner.k()
    }
    fn m(&self) -> usize {
        self.inner.m()
    }
    fn num_refs(&self) -> usize {
        self.inner.num_refs()
    }
    fn num_contigs(&self) -> usize {
        self.inner.num_contigs()
    }
    fn ref_name(&self, i: usize) -> &str {
        self.inner.ref_name(i)
    }
    fn ref_len(&self, i: usize) -> u64 {
        self.inner.ref_len(i)
    }
    fn has_ec_table(&self) -> bool {
        self.inner.has_ec_table()
    }
    fn has_poison_table(&self) -> bool {
        self.inner.has_poison_table()
    }
    fn save(&self, prefix: &str) -> Result<(), String> {
        self.inner.save(Path::new(prefix)).map_err(|e| e.to_string())
    }
}

fn make_dyn_index(index: ReferenceIndex) -> Box<dyn DynIndex> {
    let k = index.k();
    let arc = Arc::new(index);
    dispatch_on_k!(k, K => Box::new(ConcreteIndex::<K> { inner: arc }))
}

// ─── DynMappingEngine ─────────────────────────────────────────────────────────

trait DynMappingEngine: Send {
    fn map_read(&mut self, seq: &[u8]) -> MappingResultData;
    fn map_read_pair(&mut self, seq1: &[u8], seq2: &[u8]) -> MappingResultData;
    fn uses_virtual_colors(&self) -> bool;
}

struct ConcreteMappingEngine<const K: usize>
where
    Kmer<K>: KmerBits,
{
    index: Arc<ReferenceIndex>,
    strat: SkippingStrategy,
    opts: MappingOpts,
    binning: Option<BinPos>,
    cache_left: MappingCache<SketchHitInfoSimple>,
    cache_right: MappingCache<SketchHitInfoSimple>,
    cache_out: MappingCache<SketchHitInfoSimple>,
}

impl<const K: usize> ConcreteMappingEngine<K>
where
    Kmer<K>: KmerBits,
{
    fn apply_opts(cache: &mut MappingCache<SketchHitInfoSimple>, opts: &MappingOpts) {
        cache.max_hit_occ = opts.max_hit_occ;
        cache.max_read_occ = opts.max_read_occ;
        cache.max_ec_card = opts.max_ec_card;
    }

    fn extract_result(&self, cache: &MappingCache<SketchHitInfoSimple>) -> MappingResultData {
        let index = &*self.index;
        let hits: Vec<MappingHitData> = cache
            .accepted_hits
            .iter()
            .map(|h| {
                let ref_name = if (h.tid as usize) < index.num_refs() {
                    index.ref_name(h.tid as usize).to_owned()
                } else {
                    String::new()
                };
                MappingHitData {
                    tid: h.tid,
                    ref_name,
                    pos: h.pos,
                    is_fw: h.is_fw,
                    score: h.score,
                    num_kmer_hits: h.num_hits,
                    mate_pos: if h.mate_pos != INVALID_MATE_POS {
                        Some(h.mate_pos)
                    } else {
                        None
                    },
                    mate_is_fw: if h.mate_pos != INVALID_MATE_POS {
                        Some(h.mate_is_fw)
                    } else {
                        None
                    },
                    fragment_length: if h.fragment_length != INVALID_FRAG_LEN {
                        Some(h.fragment_length)
                    } else {
                        None
                    },
                    bin_id: if h.bin_id != u64::MAX {
                        Some(h.bin_id)
                    } else {
                        None
                    },
                }
            })
            .collect();

        MappingResultData {
            mapping_type: cache.map_type,
            hits,
        }
    }
}

impl<const K: usize> DynMappingEngine for ConcreteMappingEngine<K>
where
    Kmer<K>: KmerBits,
{
    fn map_read(&mut self, seq: &[u8]) -> MappingResultData {
        let index: &ReferenceIndex = &self.index;
        let mut hs = HitSearcher::new(index);
        let mut query = PiscemStreamingQuery::<K>::new(index.dict());
        let mut poison = PoisonState::new(index.poison_table());

        Self::apply_opts(&mut self.cache_out, &self.opts);

        if let Some(ref binning) = self.binning {
            map_se_fragment_atac::<K>(
                seq,
                &mut hs,
                &mut query,
                &mut self.cache_out,
                index,
                &mut poison,
                binning,
            );
        } else {
            map_se_fragment::<K>(
                seq,
                &mut hs,
                &mut query,
                &mut self.cache_out,
                index,
                &mut poison,
                self.strat,
            );
        }

        self.extract_result(&self.cache_out)
    }

    fn map_read_pair(&mut self, seq1: &[u8], seq2: &[u8]) -> MappingResultData {
        let index: &ReferenceIndex = &self.index;
        let mut hs = HitSearcher::new(index);
        let mut query = PiscemStreamingQuery::<K>::new(index.dict());
        let mut poison = PoisonState::new(index.poison_table());

        Self::apply_opts(&mut self.cache_left, &self.opts);
        Self::apply_opts(&mut self.cache_right, &self.opts);
        Self::apply_opts(&mut self.cache_out, &self.opts);

        if let Some(ref binning) = self.binning {
            map_pe_fragment_atac::<K>(
                seq1,
                seq2,
                &mut hs,
                &mut query,
                &mut self.cache_left,
                &mut self.cache_right,
                &mut self.cache_out,
                index,
                &mut poison,
                binning,
            );
        } else {
            map_pe_fragment::<K>(
                seq1,
                seq2,
                &mut hs,
                &mut query,
                &mut self.cache_left,
                &mut self.cache_right,
                &mut self.cache_out,
                index,
                &mut poison,
                self.strat,
            );
        }

        self.extract_result(&self.cache_out)
    }

    fn uses_virtual_colors(&self) -> bool {
        self.binning.is_some()
    }
}

// ─── DynStreamingQuery ────────────────────────────────────────────────────────

trait DynStreamingQuery: Send {
    fn lookup(&mut self, kmer_bytes: &[u8], read_pos: i32) -> Option<KmerHitData>;
    fn query_sequence(&mut self, seq: &[u8]) -> Vec<Option<KmerHitData>>;
    fn reset(&mut self);
    fn num_searches(&self) -> u64;
    fn num_extensions(&self) -> u64;
}

struct ConcreteStreamingQuery<const K: usize>
where
    Kmer<K>: KmerBits,
{
    query: StreamingQuery<K>,
    index: Arc<ReferenceIndex>,
    prev_query_pos: i32,
}

impl<const K: usize> ConcreteStreamingQuery<K>
where
    Kmer<K>: KmerBits,
{
    fn resolve_lookup(&self, result: &LookupResult) -> Option<KmerHitData> {
        if !result.is_found() {
            return None;
        }
        let index = &*self.index;
        let proj = index.resolve_lookup(result)?;
        let encoding = index.encoding();

        let ref_positions: Vec<RefPosData> = proj
            .ref_range()
            .iter()
            .map(|entry| {
                let rp = proj.decode_hit(entry, &encoding);
                let tid = encoding.transcript_id(entry);
                RefPosData {
                    tid,
                    pos: rp.pos,
                    is_fw: rp.is_fw,
                }
            })
            .collect();

        Some(KmerHitData {
            contig_id: proj.contig_id(),
            contig_pos: proj.contig_pos(),
            contig_len: proj.contig_len(),
            is_fw_on_contig: proj.hit_fw_on_contig(),
            ref_positions,
        })
    }
}

impl<const K: usize> DynStreamingQuery for ConcreteStreamingQuery<K>
where
    Kmer<K>: KmerBits,
{
    fn lookup(&mut self, kmer_bytes: &[u8], read_pos: i32) -> Option<KmerHitData> {
        // Auto-detect non-consecutive position and reset engine if needed
        if read_pos != self.prev_query_pos.wrapping_add(1) {
            self.query.reset();
        }
        self.prev_query_pos = read_pos;

        let result = self.query.lookup_with_dict(kmer_bytes, self.index.dict());
        self.resolve_lookup(&result)
    }

    fn query_sequence(&mut self, seq: &[u8]) -> Vec<Option<KmerHitData>> {
        let k = self.index.k();
        if seq.len() < k {
            return Vec::new();
        }
        self.query.reset();
        self.prev_query_pos = i32::MIN;

        let n = seq.len() - k + 1;
        let mut results = Vec::with_capacity(n);
        for i in 0..n {
            self.prev_query_pos = i as i32;
            let result = self.query.lookup_with_dict(&seq[i..i + k], self.index.dict());
            results.push(self.resolve_lookup(&result));
        }
        results
    }

    fn reset(&mut self) {
        self.query.reset();
        self.prev_query_pos = i32::MIN;
    }

    fn num_searches(&self) -> u64 {
        self.query.num_searches()
    }

    fn num_extensions(&self) -> u64 {
        self.query.num_extensions()
    }
}

// ─── Internal data structs ────────────────────────────────────────────────────

struct MappingResultData {
    mapping_type: MappingType,
    hits: Vec<MappingHitData>,
}

#[derive(Clone)]
struct MappingHitData {
    tid: u32,
    ref_name: String,
    pos: i32,
    is_fw: bool,
    score: f32,
    num_kmer_hits: u32,
    mate_pos: Option<i32>,
    mate_is_fw: Option<bool>,
    fragment_length: Option<i32>,
    bin_id: Option<u64>,
}

#[derive(Clone)]
struct KmerHitData {
    contig_id: u32,
    contig_pos: u32,
    contig_len: u32,
    is_fw_on_contig: bool,
    ref_positions: Vec<RefPosData>,
}

#[derive(Clone)]
struct RefPosData {
    tid: u32,
    pos: u32,
    is_fw: bool,
}

// ─── Python classes ───────────────────────────────────────────────────────────

/// A position on a reference sequence.
#[pyclass(name = "RefPos", get_all, frozen, skip_from_py_object)]
#[derive(Clone)]
struct PyRefPos {
    /// Reference ID.
    tid: u32,
    /// Position on the reference.
    pos: u32,
    /// Forward orientation on the reference.
    is_fw: bool,
}

#[pymethods]
impl PyRefPos {
    fn __repr__(&self) -> String {
        format!(
            "RefPos(tid={}, pos={}, is_fw={})",
            self.tid, self.pos, self.is_fw
        )
    }
}

/// Result of a k-mer lookup against the index.
#[pyclass(name = "KmerHit", get_all, frozen, skip_from_py_object)]
#[derive(Clone)]
struct PyKmerHit {
    /// Unitig (contig) ID.
    contig_id: u32,
    /// Position within the unitig.
    contig_pos: u32,
    /// Length of the unitig in bases.
    contig_len: u32,
    /// Whether the k-mer is in forward orientation on the unitig.
    is_fw_on_contig: bool,
    /// Decoded reference positions.
    ref_positions: Vec<PyRefPos>,
}

#[pymethods]
impl PyKmerHit {
    fn __repr__(&self) -> String {
        format!(
            "KmerHit(contig_id={}, contig_pos={}, contig_len={}, is_fw={}, num_ref_pos={})",
            self.contig_id,
            self.contig_pos,
            self.contig_len,
            self.is_fw_on_contig,
            self.ref_positions.len()
        )
    }
}

/// A single mapping hit on a reference.
#[pyclass(name = "MappingHit", get_all, frozen, skip_from_py_object)]
#[derive(Clone)]
struct PyMappingHit {
    /// Reference ID.
    tid: u32,
    /// Reference name.
    ref_name: String,
    /// Position on the reference.
    pos: i32,
    /// Forward orientation.
    is_fw: bool,
    /// Mapping score.
    score: f32,
    /// Number of supporting k-mer hits.
    num_kmer_hits: u32,
    /// Mate position (None if SE or no mate).
    mate_pos: Option<i32>,
    /// Mate orientation (None if SE or no mate).
    mate_is_fw: Option<bool>,
    /// Fragment length (None if not applicable).
    fragment_length: Option<i32>,
    /// Bin ID (only for virtual colors mode).
    bin_id: Option<u64>,
}

#[pymethods]
impl PyMappingHit {
    fn __repr__(&self) -> String {
        format!(
            "MappingHit(tid={}, ref_name='{}', pos={}, is_fw={}, score={:.1}, hits={})",
            self.tid, self.ref_name, self.pos, self.is_fw, self.score, self.num_kmer_hits
        )
    }
}

/// Result of mapping a read or read pair.
#[pyclass(name = "MappingResult", get_all, frozen)]
struct PyMappingResult {
    /// Mapping type: "unmapped", "single_mapped", "mapped_pair",
    /// "mapped_first_orphan", "mapped_second_orphan".
    mapping_type: String,
    /// Whether the read mapped.
    is_mapped: bool,
    /// List of mapping hits.
    hits: Vec<PyMappingHit>,
    /// Number of hits.
    num_hits: usize,
}

#[pymethods]
impl PyMappingResult {
    fn __repr__(&self) -> String {
        format!(
            "MappingResult(type='{}', num_hits={})",
            self.mapping_type, self.num_hits
        )
    }
}

fn mapping_type_str(mt: MappingType) -> &'static str {
    match mt {
        MappingType::Unmapped => "unmapped",
        MappingType::SingleMapped => "single_mapped",
        MappingType::MappedPair => "mapped_pair",
        MappingType::MappedFirstOrphan => "mapped_first_orphan",
        MappingType::MappedSecondOrphan => "mapped_second_orphan",
    }
}

fn to_py_mapping_result(data: MappingResultData) -> PyMappingResult {
    let is_mapped = data.mapping_type != MappingType::Unmapped;
    let hits: Vec<PyMappingHit> = data
        .hits
        .into_iter()
        .map(|h| PyMappingHit {
            tid: h.tid,
            ref_name: h.ref_name,
            pos: h.pos,
            is_fw: h.is_fw,
            score: h.score,
            num_kmer_hits: h.num_kmer_hits,
            mate_pos: h.mate_pos,
            mate_is_fw: h.mate_is_fw,
            fragment_length: h.fragment_length,
            bin_id: h.bin_id,
        })
        .collect();
    let num_hits = hits.len();
    PyMappingResult {
        mapping_type: mapping_type_str(data.mapping_type).to_owned(),
        is_mapped,
        hits,
        num_hits,
    }
}

fn to_py_kmer_hit(data: KmerHitData) -> PyKmerHit {
    let ref_positions = data
        .ref_positions
        .into_iter()
        .map(|rp| PyRefPos {
            tid: rp.tid,
            pos: rp.pos,
            is_fw: rp.is_fw,
        })
        .collect();
    PyKmerHit {
        contig_id: data.contig_id,
        contig_pos: data.contig_pos,
        contig_len: data.contig_len,
        is_fw_on_contig: data.is_fw_on_contig,
        ref_positions,
    }
}

// ─── PyMappingEngine ──────────────────────────────────────────────────────────

/// Per-read mapping engine with mutable state.
///
/// Obtain via :meth:`ReferenceIndex.mapping_engine` or
/// :meth:`ReferenceIndex.vcolor_engine`. Each Python thread should
/// create its own engine.
#[pyclass(name = "MappingEngine", unsendable)]
struct PyMappingEngine {
    engine: Box<dyn DynMappingEngine>,
}

#[pymethods]
impl PyMappingEngine {
    /// Map a single-end read.
    ///
    /// :param seq: DNA sequence as ``str`` or ``bytes``.
    /// :returns: A :class:`MappingResult`.
    fn map_read(&mut self, seq: &[u8]) -> PyMappingResult {
        to_py_mapping_result(self.engine.map_read(seq))
    }

    /// Map a paired-end read pair.
    ///
    /// :param seq1: Read 1 DNA sequence.
    /// :param seq2: Read 2 DNA sequence.
    /// :returns: A :class:`MappingResult`.
    fn map_read_pair(&mut self, seq1: &[u8], seq2: &[u8]) -> PyMappingResult {
        to_py_mapping_result(self.engine.map_read_pair(seq1, seq2))
    }

    /// Whether this engine uses virtual colors (binned mapping).
    #[getter]
    fn uses_virtual_colors(&self) -> bool {
        self.engine.uses_virtual_colors()
    }

    fn __repr__(&self) -> String {
        let mode = if self.engine.uses_virtual_colors() {
            "virtual_colors"
        } else {
            "standard"
        };
        format!("MappingEngine(mode='{mode}')")
    }
}

// ─── PyStreamingQuery ─────────────────────────────────────────────────────────

/// Low-level streaming k-mer query engine.
///
/// Maintains state across consecutive lookups for efficiency.
/// Obtain via :meth:`ReferenceIndex.streaming_query`.
#[pyclass(name = "StreamingQuery", unsendable)]
struct PyStreamingQuery {
    inner: Box<dyn DynStreamingQuery>,
}

#[pymethods]
impl PyStreamingQuery {
    /// Look up a single k-mer.
    ///
    /// :param kmer: K-mer as ``str`` or ``bytes``.
    /// :param read_pos: Position in the read (for consecutive optimization).
    /// :returns: A :class:`KmerHit` or ``None``.
    #[pyo3(signature = (kmer, read_pos=0))]
    fn lookup(&mut self, kmer: &[u8], read_pos: i32) -> Option<PyKmerHit> {
        self.inner.lookup(kmer, read_pos).map(to_py_kmer_hit)
    }

    /// Query all k-mers in a sequence.
    ///
    /// Resets the engine before processing. Returns a list where each
    /// element corresponds to one k-mer position.
    ///
    /// :param seq: DNA sequence as ``str`` or ``bytes``.
    /// :returns: ``list[KmerHit | None]``.
    fn query_sequence(&mut self, seq: &[u8]) -> Vec<Option<PyKmerHit>> {
        self.inner
            .query_sequence(seq)
            .into_iter()
            .map(|opt| opt.map(to_py_kmer_hit))
            .collect()
    }

    /// Reset the engine state. Call before processing a new sequence.
    fn reset(&mut self) {
        self.inner.reset();
    }

    /// Number of full dictionary searches performed.
    #[getter]
    fn num_searches(&self) -> u64 {
        self.inner.num_searches()
    }

    /// Number of k-mers resolved by extension.
    #[getter]
    fn num_extensions(&self) -> u64 {
        self.inner.num_extensions()
    }

    fn __repr__(&self) -> String {
        format!(
            "StreamingQuery(searches={}, extensions={})",
            self.inner.num_searches(),
            self.inner.num_extensions()
        )
    }
}

// ─── PyReferenceIndex ─────────────────────────────────────────────────────────

/// Piscem reference index for k-mer-based read mapping.
///
/// Load an existing index with :meth:`load`, or build a new one with
/// :meth:`build`. Use factory methods to create mapping engines and
/// streaming queries.
#[pyclass(name = "ReferenceIndex")]
struct PyReferenceIndex {
    inner: Box<dyn DynIndex>,
}

#[pymethods]
impl PyReferenceIndex {
    /// Load a reference index from disk.
    ///
    /// :param prefix: Index prefix (e.g., ``"path/to/gencode_v44_index"``).
    /// :param load_ec: Whether to load the equivalence class table.
    /// :param load_poison: Whether to load the poison k-mer table.
    /// :raises IOError: If the index files cannot be read.
    #[staticmethod]
    #[pyo3(signature = (prefix, *, load_ec=true, load_poison=true))]
    fn load(py: Python<'_>, prefix: &str, load_ec: bool, load_poison: bool) -> PyResult<Self> {
        let prefix_owned = prefix.to_owned();
        let index = py
            .detach(move || {
                ReferenceIndex::load(Path::new(&prefix_owned), load_ec, load_poison)
            })
            .map_err(|e| PyIOError::new_err(format!("Failed to load index '{}': {}", prefix, e)))?;
        Ok(PyReferenceIndex {
            inner: make_dyn_index(index),
        })
    }

    /// Build a new index from cuttlefish output files.
    ///
    /// :param input_prefix: Basename for cuttlefish output files (.cf_seg, .cf_seq, .json).
    /// :param output_prefix: Output prefix for index files.
    /// :param k: K-mer length (default 31).
    /// :param m: Minimizer length (default 19).
    /// :param threads: Number of threads (default 4, 0 = all cores).
    /// :param build_ec: Whether to build the equivalence class table.
    /// :param canonical: Whether to use canonical k-mers.
    /// :returns: A ready-to-use :class:`ReferenceIndex`.
    /// :raises RuntimeError: If the build fails.
    #[staticmethod]
    #[pyo3(signature = (input_prefix, output_prefix, *, k=31, m=19, threads=4, build_ec=true, canonical=true))]
    fn build(
        py: Python<'_>,
        input_prefix: &str,
        output_prefix: &str,
        k: usize,
        m: usize,
        threads: usize,
        build_ec: bool,
        canonical: bool,
    ) -> PyResult<Self> {
        use piscem_rs::index::build::BuildConfig;

        let config = BuildConfig {
            input_prefix: input_prefix.into(),
            output_prefix: output_prefix.into(),
            k,
            m,
            build_ec_table: build_ec,
            num_threads: threads,
            canonical,
            seed: 1,
            single_mphf: false,
        };
        let out_prefix = output_prefix.to_owned();
        py.detach(move || {
            piscem_rs::index::build::build_index(&config)
        })
        .map_err(|e| PyRuntimeError::new_err(format!("Build failed: {}", e)))?;

        // Load the just-built index
        let index = ReferenceIndex::load(Path::new(&out_prefix), true, false)
            .map_err(|e| {
                PyRuntimeError::new_err(format!("Failed to load built index: {}", e))
            })?;
        Ok(PyReferenceIndex {
            inner: make_dyn_index(index),
        })
    }

    /// Save the index to disk.
    ///
    /// :param prefix: Output prefix.
    fn save(&self, prefix: &str) -> PyResult<()> {
        self.inner
            .save(prefix)
            .map_err(|e| PyIOError::new_err(format!("Failed to save index: {}", e)))
    }

    /// K-mer size.
    #[getter]
    fn k(&self) -> usize {
        self.inner.k()
    }

    /// Minimizer length.
    #[getter]
    fn m(&self) -> usize {
        self.inner.m()
    }

    /// Number of reference sequences.
    #[getter]
    fn num_refs(&self) -> usize {
        self.inner.num_refs()
    }

    /// Number of unitigs (contigs).
    #[getter]
    fn num_contigs(&self) -> usize {
        self.inner.num_contigs()
    }

    /// Whether an equivalence class table is loaded.
    #[getter]
    fn has_ec_table(&self) -> bool {
        self.inner.has_ec_table()
    }

    /// Whether a poison k-mer table is loaded.
    #[getter]
    fn has_poison_table(&self) -> bool {
        self.inner.has_poison_table()
    }

    /// Get the name of reference ``i``.
    fn ref_name(&self, i: usize) -> PyResult<String> {
        if i >= self.inner.num_refs() {
            return Err(PyValueError::new_err(format!(
                "Reference index {} out of range (num_refs={})",
                i,
                self.inner.num_refs()
            )));
        }
        Ok(self.inner.ref_name(i).to_owned())
    }

    /// Get the length of reference ``i``.
    fn ref_len(&self, i: usize) -> PyResult<u64> {
        if i >= self.inner.num_refs() {
            return Err(PyValueError::new_err(format!(
                "Reference index {} out of range (num_refs={})",
                i,
                self.inner.num_refs()
            )));
        }
        Ok(self.inner.ref_len(i))
    }

    /// Get all reference names.
    fn ref_names(&self) -> Vec<String> {
        (0..self.inner.num_refs())
            .map(|i| self.inner.ref_name(i).to_owned())
            .collect()
    }

    /// Get all reference lengths.
    fn ref_lengths(&self) -> Vec<u64> {
        (0..self.inner.num_refs())
            .map(|i| self.inner.ref_len(i))
            .collect()
    }

    /// Create a mapping engine for standard read mapping.
    ///
    /// :param strategy: ``"permissive"`` (default) or ``"strict"``.
    /// :param max_hit_occ: Maximum k-mer occurrence threshold (default 256).
    /// :param max_read_occ: Maximum mappings per read before discarding (default 2500).
    /// :returns: A :class:`MappingEngine`.
    #[pyo3(signature = (*, strategy="permissive", max_hit_occ=256, max_read_occ=2500))]
    fn mapping_engine(
        &self,
        strategy: &str,
        max_hit_occ: usize,
        max_read_occ: usize,
    ) -> PyResult<PyMappingEngine> {
        let strat = match strategy {
            "permissive" => SkippingStrategy::Permissive,
            "strict" => SkippingStrategy::Strict,
            _ => {
                return Err(PyValueError::new_err(format!(
                    "Unknown strategy '{}': use 'permissive' or 'strict'",
                    strategy
                )));
            }
        };
        let opts = MappingOpts {
            max_hit_occ,
            max_read_occ,
            ..MappingOpts::default()
        };
        Ok(PyMappingEngine {
            engine: self.inner.make_mapping_engine(strat, opts),
        })
    }

    /// Create a mapping engine using virtual colors (binned mapping).
    ///
    /// Used for scATAC-seq and similar binned genomic mapping.
    ///
    /// :param bin_size: Size of each bin in bases (default 2000).
    /// :param overlap: Overlap between adjacent bins (default 400).
    /// :param thr: Hit threshold fraction (default 0.7).
    /// :param max_hit_occ: Maximum k-mer occurrence threshold (default 256).
    /// :returns: A :class:`MappingEngine`.
    #[pyo3(signature = (*, bin_size=2000, overlap=400, thr=0.7, max_hit_occ=256))]
    fn vcolor_engine(
        &self,
        bin_size: u64,
        overlap: u64,
        thr: f32,
        max_hit_occ: usize,
    ) -> PyMappingEngine {
        let opts = MappingOpts {
            max_hit_occ,
            ..MappingOpts::default()
        };
        PyMappingEngine {
            engine: self.inner.make_vcolor_engine(opts, bin_size, overlap, thr),
        }
    }

    /// Create a low-level streaming k-mer query engine.
    ///
    /// :returns: A :class:`StreamingQuery`.
    fn streaming_query(&self) -> PyStreamingQuery {
        PyStreamingQuery {
            inner: self.inner.make_streaming_query(),
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "ReferenceIndex(k={}, m={}, num_refs={}, num_contigs={}, ec={}, poison={})",
            self.inner.k(),
            self.inner.m(),
            self.inner.num_refs(),
            self.inner.num_contigs(),
            self.inner.has_ec_table(),
            self.inner.has_poison_table(),
        )
    }
}

// ─── Module ───────────────────────────────────────────────────────────────────

/// Python bindings for the piscem k-mer mapping engine.
#[pymodule]
fn piscem(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyReferenceIndex>()?;
    m.add_class::<PyMappingEngine>()?;
    m.add_class::<PyMappingResult>()?;
    m.add_class::<PyMappingHit>()?;
    m.add_class::<PyStreamingQuery>()?;
    m.add_class::<PyKmerHit>()?;
    m.add_class::<PyRefPos>()?;
    Ok(())
}
