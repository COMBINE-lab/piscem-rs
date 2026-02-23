"""Tests for the piscem Python bindings.

Requires a pre-built index. Set PISCEM_TEST_INDEX to the index prefix, e.g.:
    PISCEM_TEST_INDEX=/path/to/gencode_pc_v44_index pytest tests/
"""

import os
import gzip
import pytest
import piscem

INDEX_PREFIX = os.environ.get(
    "PISCEM_TEST_INDEX",
    os.path.join(
        os.path.dirname(__file__),
        "..", "..", "..",
        "test_data", "gencode_pc_v44_index_rust", "gencode_pc_v44_index",
    ),
)

READS_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "..", "test_data")


def _have_index():
    return os.path.exists(INDEX_PREFIX + ".ssi")


def _have_reads():
    return os.path.exists(os.path.join(READS_DIR, "sim_1M_1.fq.gz"))


def _read_first_pair():
    """Read the first read pair from the test FASTQ files."""
    with gzip.open(os.path.join(READS_DIR, "sim_1M_1.fq.gz"), "rt") as f:
        f.readline()
        seq1 = f.readline().strip()
    with gzip.open(os.path.join(READS_DIR, "sim_1M_2.fq.gz"), "rt") as f:
        f.readline()
        seq2 = f.readline().strip()
    return seq1, seq2


@pytest.fixture(scope="module")
def index():
    if not _have_index():
        pytest.skip("Test index not found; set PISCEM_TEST_INDEX")
    return piscem.ReferenceIndex.load(INDEX_PREFIX)


# ── Index loading and metadata ───────────────────────────────────────────────

class TestIndexMetadata:
    def test_load(self, index):
        assert index.k == 31
        assert index.m == 19
        assert index.num_refs > 0
        assert index.num_contigs > 0

    def test_repr(self, index):
        r = repr(index)
        assert "ReferenceIndex" in r
        assert "k=31" in r

    def test_ref_name(self, index):
        name = index.ref_name(0)
        assert isinstance(name, str)
        assert len(name) > 0

    def test_ref_len(self, index):
        length = index.ref_len(0)
        assert isinstance(length, int)
        assert length > 0

    def test_ref_names_lengths(self, index):
        names = index.ref_names()
        lengths = index.ref_lengths()
        assert len(names) == index.num_refs
        assert len(lengths) == index.num_refs
        assert names[0] == index.ref_name(0)
        assert lengths[0] == index.ref_len(0)

    def test_ref_name_bounds(self, index):
        with pytest.raises(ValueError):
            index.ref_name(index.num_refs)

    def test_ref_len_bounds(self, index):
        with pytest.raises(ValueError):
            index.ref_len(index.num_refs)

    def test_has_ec_table(self, index):
        assert isinstance(index.has_ec_table, bool)

    def test_has_poison_table(self, index):
        assert isinstance(index.has_poison_table, bool)


# ── Mapping engine ───────────────────────────────────────────────────────────

class TestMappingEngine:
    def test_create_engine(self, index):
        eng = index.mapping_engine()
        assert repr(eng) == "MappingEngine(mode='standard')"
        assert not eng.uses_virtual_colors

    def test_create_strict_engine(self, index):
        eng = index.mapping_engine(strategy="strict")
        assert not eng.uses_virtual_colors

    def test_invalid_strategy(self, index):
        with pytest.raises(ValueError):
            index.mapping_engine(strategy="bogus")

    def test_map_unmapped_read(self, index):
        eng = index.mapping_engine()
        result = eng.map_read(b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
        assert result.mapping_type == "unmapped"
        assert not result.is_mapped
        assert result.num_hits == 0
        assert result.hits == []

    @pytest.mark.skipif(not _have_reads(), reason="Test reads not found")
    def test_map_real_se(self, index):
        seq1, _ = _read_first_pair()
        eng = index.mapping_engine()
        result = eng.map_read(seq1.encode())
        assert result.is_mapped
        assert result.mapping_type == "single_mapped"
        assert result.num_hits > 0
        hit = result.hits[0]
        assert isinstance(hit.tid, int)
        assert isinstance(hit.ref_name, str)
        assert isinstance(hit.pos, int)
        assert isinstance(hit.is_fw, bool)
        assert isinstance(hit.score, float)
        assert isinstance(hit.num_kmer_hits, int)

    @pytest.mark.skipif(not _have_reads(), reason="Test reads not found")
    def test_map_real_pe(self, index):
        seq1, seq2 = _read_first_pair()
        eng = index.mapping_engine()
        result = eng.map_read_pair(seq1.encode(), seq2.encode())
        assert result.is_mapped
        assert result.mapping_type == "mapped_pair"
        hit = result.hits[0]
        assert hit.fragment_length is not None
        assert hit.mate_pos is not None

    def test_map_accepts_bytes(self, index):
        eng = index.mapping_engine()
        result = eng.map_read(b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
        assert result.mapping_type == "unmapped"

    def test_create_struct_constraints_engine(self, index):
        eng = index.mapping_engine(struct_constraints=True)
        assert not eng.uses_virtual_colors

    @pytest.mark.skipif(not _have_reads(), reason="Test reads not found")
    def test_map_struct_constraints_se(self, index):
        seq1, _ = _read_first_pair()
        eng = index.mapping_engine(struct_constraints=True)
        result = eng.map_read(seq1.encode())
        assert result.is_mapped
        assert result.mapping_type == "single_mapped"

    @pytest.mark.skipif(not _have_reads(), reason="Test reads not found")
    def test_map_struct_constraints_pe(self, index):
        seq1, seq2 = _read_first_pair()
        eng = index.mapping_engine(struct_constraints=True)
        result = eng.map_read_pair(seq1.encode(), seq2.encode())
        assert result.is_mapped


# ── Virtual colors engine ────────────────────────────────────────────────────

class TestVColorEngine:
    def test_create_vcolor_engine(self, index):
        eng = index.vcolor_engine()
        assert eng.uses_virtual_colors
        assert "virtual_colors" in repr(eng)

    @pytest.mark.skipif(not _have_reads(), reason="Test reads not found")
    def test_vcolor_map(self, index):
        seq1, seq2 = _read_first_pair()
        eng = index.vcolor_engine(bin_size=2000, overlap=400, thr=0.7)
        result = eng.map_read_pair(seq1.encode(), seq2.encode())
        # May or may not map depending on binning params
        assert isinstance(result.mapping_type, str)


# ── Streaming query ──────────────────────────────────────────────────────────

class TestStreamingQuery:
    def test_create(self, index):
        sq = index.streaming_query()
        assert sq.num_searches == 0
        assert sq.num_extensions == 0

    def test_query_short_seq(self, index):
        sq = index.streaming_query()
        # Shorter than k → empty list
        hits = sq.query_sequence(b"ACGT")
        assert hits == []

    @pytest.mark.skipif(not _have_reads(), reason="Test reads not found")
    def test_query_real_sequence(self, index):
        seq1, _ = _read_first_pair()
        sq = index.streaming_query()
        hits = sq.query_sequence(seq1.encode())
        assert len(hits) == len(seq1) - index.k + 1
        found = [h for h in hits if h is not None]
        assert len(found) > 0

        hit = found[0]
        assert isinstance(hit.contig_id, int)
        assert isinstance(hit.contig_pos, int)
        assert isinstance(hit.contig_len, int)
        assert isinstance(hit.is_fw_on_contig, bool)
        assert len(hit.ref_positions) > 0

        rp = hit.ref_positions[0]
        assert isinstance(rp.tid, int)
        assert isinstance(rp.pos, int)
        assert isinstance(rp.is_fw, bool)

    def test_reset(self, index):
        sq = index.streaming_query()
        sq.reset()
        assert sq.num_searches == 0

    def test_lookup(self, index):
        sq = index.streaming_query()
        # Random k-mer likely not found
        result = sq.lookup(b"A" * index.k)
        # Could be None or a hit; just test it doesn't crash
        assert result is None or isinstance(result, piscem.KmerHit)
