"""Type stubs for the piscem Python bindings."""

from __future__ import annotations

class RefPos:
    """A position on a reference sequence."""

    tid: int
    pos: int
    is_fw: bool

class KmerHit:
    """Result of a k-mer lookup against the index."""

    contig_id: int
    contig_pos: int
    contig_len: int
    is_fw_on_contig: bool
    ref_positions: list[RefPos]

class MappingHit:
    """A single mapping hit on a reference."""

    tid: int
    ref_name: str
    pos: int
    is_fw: bool
    score: float
    num_kmer_hits: int
    mate_pos: int | None
    mate_is_fw: bool | None
    fragment_length: int | None
    bin_id: int | None

class MappingResult:
    """Result of mapping a read or read pair."""

    mapping_type: str
    is_mapped: bool
    hits: list[MappingHit]
    num_hits: int

class MappingEngine:
    """Per-read mapping engine with mutable state.

    Obtain via :meth:`ReferenceIndex.mapping_engine` or
    :meth:`ReferenceIndex.vcolor_engine`.
    """

    uses_virtual_colors: bool

    def map_read(self, seq: bytes) -> MappingResult: ...
    def map_read_pair(self, seq1: bytes, seq2: bytes) -> MappingResult: ...

class StreamingQuery:
    """Low-level streaming k-mer query engine.

    Obtain via :meth:`ReferenceIndex.streaming_query`.
    """

    num_searches: int
    num_extensions: int

    def lookup(self, kmer: bytes, read_pos: int = 0) -> KmerHit | None: ...
    def query_sequence(self, seq: bytes) -> list[KmerHit | None]: ...
    def reset(self) -> None: ...

class ReferenceIndex:
    """Piscem reference index for k-mer-based read mapping."""

    k: int
    m: int
    num_refs: int
    num_contigs: int
    has_ec_table: bool
    has_poison_table: bool

    @staticmethod
    def load(
        prefix: str,
        *,
        load_ec: bool = True,
        load_poison: bool = True,
    ) -> ReferenceIndex: ...
    @staticmethod
    def build(
        input_prefix: str,
        output_prefix: str,
        *,
        k: int = 31,
        m: int = 19,
        threads: int = 4,
        build_ec: bool = True,
        canonical: bool = True,
        decoys: list[str] | None = None,
    ) -> ReferenceIndex: ...
    def save(self, prefix: str) -> None: ...
    def ref_name(self, i: int) -> str: ...
    def ref_len(self, i: int) -> int: ...
    def ref_names(self) -> list[str]: ...
    def ref_lengths(self) -> list[int]: ...
    def mapping_engine(
        self,
        *,
        strategy: str = "permissive",
        max_hit_occ: int = 256,
        max_read_occ: int = 2500,
        struct_constraints: bool = False,
    ) -> MappingEngine: ...
    def vcolor_engine(
        self,
        *,
        bin_size: int = 2000,
        overlap: int = 400,
        thr: float = 0.7,
        max_hit_occ: int = 256,
    ) -> MappingEngine: ...
    def streaming_query(self) -> StreamingQuery: ...
