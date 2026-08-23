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

def dedupe_reads(
    inputs: List[str],
    min_length: int,
    batch_size: int,
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
def filter_nt(candidates: List[Tuple[str, str]], failed: Set[str]) -> List[Tuple[str, str]]: ...

# Column operations
def del_cols(sequence: str, x_positions: Set[int], is_nt: bool) -> str: ...
def join_by_tripled_index(string: str, positions_to_keep: List[int]) -> str: ...
def join_with_exclusions(string: str, column_cull: Set[int]) -> str: ...
def join_triplets_with_exclusions(
    string: str,
    exclusion1: Set[int],
    exclusion2: Set[int],
) -> str: ...
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
) -> List[Tuple[str, str]]: ...
def exon_dp(folder: str, sub_dir: str, taxa_path: str) -> Any: ...
# flank_inputs, gap_inputs and gff_nodes are extracted by attribute name
# (FlankInput / GapInput / GffEntryInput), so any object with the right
# attributes is accepted.
def exonfinder_process_gene(
    gene_key: str,
    aa_file: str,
    refs_aa: List[Tuple[str, str]],
    natives_aa: List[Tuple[str, str]],
    flank_inputs: List[Any],
    gap_inputs: List[Any],
    nt_seqs: Dict[str, str],
    gff_nodes: Dict[str, Any],
    clusters_for_gene: Dict[str, List[str]],
    taxa: Optional[str],
    tmpdir: Optional[str],
    debug_level: int,
    skip_exon_split: bool,
    disable_column_cull: bool,
    log_dir: Optional[str],
    disable_score_filter: bool,
) -> Dict[str, Any]: ...
