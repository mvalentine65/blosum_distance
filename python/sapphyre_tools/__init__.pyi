from typing import Any, Dict, List, Optional, Set, Tuple

class NtBatchScanner:
    def __init__(self, nodes: List[int]) -> None: ...
    def __len__(self) -> int: ...
    def scan_into(self, batch: bytes, out: Dict[int, bytes]) -> int: ...

class PreparedReads:
    def record_count(self) -> int: ...
    def total_dupes(self) -> int: ...
    def batch_count(self) -> int: ...
    def batch(self, index: int) -> bytes: ...
    def packed_dupes(self) -> bytes: ...
    def reads_in(self) -> int: ...
    def reads_passed(self) -> int: ...
    def adapter_trimmed_reads(self) -> int: ...
    def adapter_trimmed_bases(self) -> int: ...
    def filtered_counts(self) -> Tuple[int, int, int, int]: ...
    def detected_adapters(self) -> Tuple[Optional[str], Optional[str]]: ...
    def base_counts(self) -> Tuple[int, int]: ...
    def polyx_trimmed(self) -> Tuple[int, int]: ...
    def length_histogram(self) -> List[int]: ...
    def repeats_removed(self) -> Tuple[int, int, int]: ...
    def repeat_units(self) -> Tuple[List[Tuple[str, int]], int]: ...
    def insert_size_histogram(self) -> List[int]: ...

class CullTables:
    def __init__(
        self,
        character_at_each_pos: Dict[int, Any],
        blosum_at_each_pos: Dict[int, Any],
        all_dashes_by_index: Dict[int, Any],
        gap_present_threshold: Dict[int, Any],
    ) -> None: ...
    def __len__(self) -> int: ...
    def cull_many(
        self,
        sequences: List[str],
        offset: int,
        amt_matches: int,
        mismatches: int,
        blosum_max_percent: float,
    ) -> List[Tuple[Optional[int], Optional[int], bool]]: ...

def dedupe_reads(
    inputs: List[str],
    min_length: int,
    batch_size: int,
    trim: bool,
    repeat_filter: bool,
) -> PreparedReads: ...

# Sequence distance / scoring
def blosum62_distance(one: str, two: str) -> float: ...
def blosum62_candidate_to_reference(candidate: str, reference: str) -> float: ...
def constrained_distance(consensus: str, candidate: str) -> int: ...
def consensus_distance(
    consensus: str,
    candidate: str,
    min_length: int,
    min_overlap: int,
) -> Tuple[int, int]: ...

# Consensus
def dumb_consensus(sequences: List[str], threshold: float, min_depth: int) -> str: ...
def dumb_consensus_dupe(
    sequences: List[Tuple[str, int]],
    threshold: float,
    min_depth: int,
) -> str: ...
def convert_consensus(sequences: List[str], consensus: str) -> str: ...

# Sequence utilities
def bio_revcomp(sequence: str) -> str: ...
def translate(sequence: str, table: Optional[int] = ...) -> str: ...
def find_index_pair(sequence: str, gap: str) -> Tuple[int, int]: ...
def is_same_kmer(one: str, two: str) -> bool: ...
def get_overlap(
    start1: int,
    end1: int,
    start2: int,
    end2: int,
    min_overlap: int,
) -> Optional[Tuple[int, int]]: ...
def is_low_complexity_nt(
    seq: str,
    window_size: int,
    step: int,
    dinuc_entropy_threshold: float,
    trinuc_entropy_threshold: float,
    min_segment_length: int,
    min_failing_fraction: float,
    min_failing_bp: int,
    min_total_windows: int,
) -> bool: ...

# Column operations
def join_with_exclusions(string: str, column_cull: Set[int]) -> str: ...
def join_triplets_with_exclusions_many(
    sequences: List[str],
    exclusions: List[Set[int]],
    shared_exclusion: Set[int],
) -> List[str]: ...
def delete_empty_columns_pairs(
    records: List[Tuple[str, str]],
) -> List[Tuple[str, str]]: ...
def cull_columns(
    records: List[Tuple[str, str]],
    ref_suffix: str,
    max_allowed_gaps_in_ref: float,
    min_ref_supported_gap: int,
    anchor_run: int,
    tail_max_data: int,
    include_edge: bool,
    edge_ref_data_fraction: float,
    min_seq_len: int,
    gap_cull_threshold: float,
    nt_seqs: Optional[Dict[str, str]],
    debug: int,
    log_dir: Optional[str],
) -> Tuple[
    List[Tuple[str, str]],
    Dict[str, Set[int]],
    Dict[str, Tuple[int, int]],
    Dict[str, List[Tuple[int, int]]],
    Dict[str, List[Tuple[int, int]]],
]: ...
def apply_gff_culls(
    source_gff_path: str,
    output_gff_path: str,
    culls: Dict[Tuple[str, str, Optional[str]], Optional[Tuple[int, int]]],
    intron_splits: Optional[
        Dict[Tuple[str, str, Optional[str]], List[Tuple[int, int]]]
    ] = ...,
) -> Tuple[bool, Set[str]]: ...

# Alignment / exon finding
def hmm_align(
    candidates: List[Tuple[str, str]],
    references: List[Tuple[str, str]],
    tmpdir: Optional[str] = ...,
    gene_name: Optional[str] = ...,
    taxa: Optional[str] = ...,
    cached_hmm: Optional[str] = ...,
    cached_template: Optional[str] = ...,
) -> List[Tuple[str, str]]: ...

# Codon alignment
def pn2codon(
    _file_steem: str,
    aa_path: str,
    nt_path: str,
    table_num: int,
    seqs: Dict[str, Tuple[Tuple[str, str], Tuple[int, str, str]]],
) -> str: ...

# Exon recovery
def exonfill_model_cols(models: str, threads: int = ...) -> Dict[str, List[int]]: ...
def exonfill_run(
    models: str,
    chains: str,
    genome: Any,
    prefix: str,
    opts: Optional[Dict[str, float]] = ...,
    stages: str = ...,
) -> None: ...
