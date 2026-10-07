//! What joins a row at its ends, decided in one place: whether a recovered ORF
//! piece is kept, and how far a row's first and last exon may reach for a start
//! or a stop codon. Each rule is handed all the evidence there is about the
//! piece or the end, whether it reads it yet or not.

use super::orf::OrfOpts;

/// A recovered ORF piece as it is judged. `cover` is evidence no rule reads yet.
#[allow(dead_code)]
#[derive(Clone, Copy)]
pub struct Piece {
    /// Forward bits against the window's model nodes
    pub bits: f32,
    /// bits less the bits of the piece reversed
    pub margin: f32,
    /// bits charged for the ORFs tried in the window
    pub charge: f32,
    /// residues on model nodes
    pub nmatch: i32,
    /// codons aligned
    pub aa: i64,
    pub repeat: bool,
    /// nt between the piece and the row's exon beside the window (the nearer one in a gap)
    pub dist: i64,
    /// model nodes between the piece and that exon (the fewer in a gap)
    pub skip: i64,
    /// found past a row's end, not between two of its exons
    pub flank: bool,
    /// share of the model's nodes the row's exons cover
    pub cover: f64,
    /// found before the row's first exon
    pub lead: bool,
    /// that exon reads back to a start codon of its own (lead_start_reach)
    pub start: bool,
    /// in a gap: the two exons beside it join on their own (fill_gap)
    pub joined: bool,
}

/// nt within which a piece has no room for an intron before the row's exon
const CLOSE_NT: i64 = 30;
/// nt beyond which a piece is far from the row
const FAR_NT: i64 = 100;
/// codons under which a piece is short
const SHORT_AA: i64 = 20;
/// model nodes a flank piece may lie from its exon before it needs FAR_NODES_BITS
const FAR_NODES: i64 = 50;
const FAR_NODES_BITS: f32 = 10.0;

/// Bits a piece must score for where it lies: the less its place supports it,
/// the more. A piece 31 to 100 nt from the row's exon, past the row's end and
/// of 20 codons or more, needs none. One with no room for an intron needs 2,
/// one far from the row 3; between two exons 2 more, short 3 more. A flank
/// piece over FAR_NODES model nodes from its exon needs FAR_NODES_BITS.
pub fn bits_needed(p: &Piece) -> f32 {
    let mut need = 0.0;
    if p.dist <= CLOSE_NT { need += 2.0; }
    if p.dist > FAR_NT { need += 3.0; }
    if !p.flank { need += 2.0; }
    if p.aa < SHORT_AA { need += 3.0; }
    if p.flank && p.skip > FAR_NODES { need = f32::max(need, FAR_NODES_BITS); }
    need
}

/// nt from which a piece is far before a first exon that has its own start codon
const LEAD_FAR_NT: i64 = 1000;
/// codons under which such a piece is short
const LEAD_SHORT_AA: i64 = 30;

/// A short piece far before a first exon that reads back to a start codon of
/// its own owes its intron: the gene may as well start at that codon, and the
/// model's first nodes match short hydrophobic ORFs easily. The search keeps it;
/// it is judged once its junction has sites (pays_intron).
pub fn owes_intron(p: &Piece) -> bool {
    p.lead && p.start && p.dist >= LEAD_FAR_NT && p.aa < LEAD_SHORT_AA
}

/// A piece that owes its intron stays when its bits outweigh what the spliced
/// search charges for that intron (`intron`, half-bits: splice::intron_charge).
pub fn pays_intron(bits: f32, intron: f64) -> bool {
    bits as f64 + intron / 2.0 > 0.0
}

/// A piece in a gap whose two exons already join on their own, with no model
/// node between them left over, is refused: they cover what it would, and in
/// the row it forces an intron on each side of itself.
pub fn gap_covered(p: &Piece) -> bool {
    !p.flank && p.joined
}

/// A piece is kept when it is no repeat, beats its reversal by `revthr`, has
/// `minm` residues on nodes, either outscores the charge by `thr` or lies
/// within `near` nt of the row, scores what its place asks (bits_needed), and
/// is no covered gap piece.
pub fn keep_piece(p: &Piece, o: &OrfOpts) -> bool {
    !p.repeat && p.margin >= o.revthr && p.nmatch >= o.minm && (p.bits - p.charge >= o.thr || p.dist <= o.near)
        && p.bits >= bits_needed(p) && !gap_covered(p)
}

/// nt the pieces of a run may lie apart, and the first from the row's exon
const RUN_NT: i64 = 1000;
/// model nodes two pieces of a run may share
const RUN_SHARE: i64 = 3;
/// bits the pieces of a run sum to
const RUN_BITS: f32 = 10.0;

/// A piece may be one of a run when it is no repeat, beats its reversal by
/// `revthr` and has `minm` residues on nodes; its own bits are not asked.
pub fn in_run(p: &Piece, o: &OrfOpts) -> bool {
    !p.repeat && p.margin >= o.revthr && p.nmatch >= o.minm && p.bits > 0.0
}

/// Pieces past a row's end, none kept alone, that follow one another in the
/// genome and in the model (run_next), are judged as one piece: `n` of them,
/// two or more, summing to RUN_BITS and outscoring the charge by `thr`. The run
/// starts within RUN_NT and FAR_NODES of the row's exon (run_start). Exons of
/// one gene each too short or too far in the model to stand alone.
pub fn keep_run(bits: f32, n: usize, charge: f32, o: &OrfOpts) -> bool {
    n >= 2 && bits >= RUN_BITS && bits - charge >= o.thr
}

/// Whether a run may start at this piece, the one next to the row's exon.
pub fn run_start(p: &Piece) -> bool {
    p.dist <= RUN_NT && p.skip <= FAR_NODES
}

/// Whether a piece follows the one before it in a run: `gap` nt past it with
/// room for an intron, `skip` model nodes on (negative: nodes shared).
pub fn run_next(gap: i64, skip: i64, min_intron: i64) -> bool {
    gap >= min_intron && gap <= RUN_NT && skip >= -RUN_SHARE && skip <= FAR_NODES
}

/// A piece kept as one of a run must meet both neighbours across a junction
/// with splice sites and no frameshift or stop (`clean`), and the run must still
/// hold every piece it had (`whole`); else the run goes.
pub fn keep_run_piece(clean: bool, whole: bool) -> bool {
    clean && whole
}

/// A piece that meets a neighbour only across a strict frameshift needs `floor`
/// bits: under that it is intron read as exon.
pub fn keep_at_frameshift(bits: f32, floor: f32) -> bool {
    bits >= floor
}

/// An internal recovered exon must fit its chain at least `min_rev` bits better
/// than reversed. NaN: untested, kept.
pub fn keep_against_reversal(rev: f64, min_rev: f64) -> bool {
    !(rev < min_rev)
}

/// A row's end as its extension is judged. `exons` and `lone` are evidence no
/// rule reads yet.
#[allow(dead_code)]
pub struct End {
    /// model nodes between the end exon and that end of the model
    pub short: i64,
    /// model nodes the row's exons cover
    pub covered: usize,
    /// the model's length
    pub m: usize,
    /// exons in the row
    pub exons: usize,
    /// a hit in no row
    pub lone: bool,
}

const START_NODES: i64 = 30;
const START_NT: i64 = 150;
const START_COVER: f64 = 0.7;
/// codons a start codon may lie beyond the model nodes the first exon lacks
const START_SLACK: i64 = 20;
const STOP_NODES: i64 = 20;
const STOP_NT: i64 = 90;
const STOP_REACH_NODES: i64 = 100;
/// nt past a row's last exon read for a stop codon
pub const STOP_SCAN_NT: i64 = 1500;

/// nt upstream of a row's first exon to look for its start codon (the farthest
/// in-frame ATG before a stop), when the row earns the search: the exon starts
/// within START_NODES of the model's first node and the row's exons cover
/// START_COVER of the model. A lone hit or a thin row is too often no gene's
/// first exon. Genes run on before the model's first node. START_NT covers the
/// reach lead_start_reach looks in.
pub fn start_reach(e: &End) -> Option<i64> {
    (e.short < START_NODES && e.covered as f64 >= START_COVER * e.m as f64).then_some(START_NT)
}

/// nt upstream of a first exon, `short` model nodes from the model's start, in
/// which a start codon of its own counts against a far lead piece: the nodes it
/// lacks and START_SLACK codons. None when the exon is not near the model's start.
pub fn lead_start_reach(short: i64) -> Option<i64> {
    (short < START_NODES).then_some(3 * (short + START_SLACK))
}

/// Whether the first in-frame stop, `dist` nt past a row's last exon, ends the
/// row. Within STOP_NODES of the model's last node a stop within STOP_NT does.
/// Further from the model's end, up to STOP_REACH_NODES, any stop read does
/// unless a donor lies on the way (`donor`): a further exon. Beyond that the
/// exon is rarely the gene's last and is left as it is.
pub fn take_stop(e: &End, dist: i64, donor: impl FnOnce() -> bool) -> bool {
    (e.short <= STOP_NODES && dist <= STOP_NT) || (e.short <= STOP_REACH_NODES && !donor())
}
