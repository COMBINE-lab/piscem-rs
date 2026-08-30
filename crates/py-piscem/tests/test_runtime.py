"""Runtime tests for the piscem python bindings.

Run with: maturin develop --release && pytest crates/py-piscem/tests

Data-gated like the workspace's RAD parity tests: set
  PISCEM_PY_INDEX     — piscem index prefix (e.g. .../index)
  PISCEM_PY_TXP_FASTA — the transcripts FASTA the index was built from
  PISCEM_PY_READS1/2  — paired FASTQ files whose reads map to those transcripts
  PISCEM_PY_CF_PREFIX — cuttlefish output prefix (.cf_seg/.cf_seq) [optional,
                        enables the build() test]
otherwise the tests are skipped. The deprecation-warning test for the
`canonical=` kwarg only needs PISCEM_PY_CF_PREFIX.
"""
import os
import warnings

import pytest

piscem = pytest.importorskip("piscem")

IDX = os.environ.get("PISCEM_PY_INDEX")
FASTA = os.environ.get("PISCEM_PY_TXP_FASTA")
R1 = os.environ.get("PISCEM_PY_READS1")
R2 = os.environ.get("PISCEM_PY_READS2")
CF = os.environ.get("PISCEM_PY_CF_PREFIX")

needs_index = pytest.mark.skipif(
    not (IDX and FASTA), reason="PISCEM_PY_INDEX / PISCEM_PY_TXP_FASTA not set")


def read_fasta(path):
    recs, name, cur = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                recs[name] = "".join(cur)
            name, cur = line[1:].split()[0], []
        else:
            cur.append(line.strip().upper())
    if name:
        recs[name] = "".join(cur)
    return recs


def read_fastq(path, n):
    reads, lines = [], open(path).read().splitlines()
    for i in range(0, min(4 * n, len(lines)), 4):
        reads.append(lines[i + 1].upper())
    return reads


@pytest.fixture(scope="module")
def idx():
    return piscem.ReferenceIndex.load(IDX)


@pytest.fixture(scope="module")
def txps():
    return read_fasta(FASTA)


@needs_index
def test_shape(idx, txps):
    assert idx.num_refs == len(txps)
    assert set(idx.ref_names()) == set(txps)
    names = idx.ref_names()
    for i in range(idx.num_refs):
        assert idx.ref_len(i) <= len(txps[names[i]])  # poly-A clipping may shorten


@needs_index
def test_streaming_query_finds_indexed_kmers(idx, txps):
    k = idx.k
    seq = txps[idx.ref_names()[0]]
    sq = idx.streaming_query()
    n = min(300, len(seq) - k + 1)
    for i in range(n):
        assert sq.lookup(seq[i : i + k].encode()) is not None, f"miss at {i}"


@pytest.mark.skipif(not (IDX and R1 and R2), reason="reads not set")
def test_mapping_engine_maps_reads(idx):
    eng = idx.mapping_engine()
    r1 = read_fastq(R1, 200)
    r2 = read_fastq(R2, 200)
    mapped = 0
    for a, b in zip(r1, r2):
        res = eng.map_read_pair(a.encode(), b.encode())
        if res.is_mapped:
            mapped += 1
            for h in res.hits:
                assert 0 <= h.tid < idx.num_refs
    assert mapped > 0.9 * len(r1)


@pytest.mark.skipif(not CF, reason="PISCEM_PY_CF_PREFIX not set")
def test_build_canonical_kwarg_deprecation(tmp_path):
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        built = piscem.ReferenceIndex.build(
            CF, str(tmp_path / "pyidx"),
            k=31, m=19, threads=2, build_ec=False, canonical=True)
    dep = [x for x in w if issubclass(x.category, DeprecationWarning)]
    assert len(dep) == 1
    assert "always canonical" in str(dep[0].message)
    assert built.k == 31
    # omitted kwarg is silent
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        piscem.ReferenceIndex.build(
            CF, str(tmp_path / "pyidx2"),
            k=31, m=19, threads=2, build_ec=False)
    assert not [x for x in w if issubclass(x.category, DeprecationWarning)]
