//! What joins a row at its ends, decided in one place: whether a recovered ORF
//! piece is kept, and how far a row's first and last exon may reach for a start
//! or a stop codon. Each rule is handed all the evidence there is about the
//! piece or the end, whether it reads it yet or not.

use super::orf::OrfOpts;

/// A recovered ORF piece as it is judged. `cover` is evidence no rule reads yet.
#[allow(dead_code)]
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

/// A piece is kept when it is no repeat, beats its reversal by `revthr`, has
/// `minm` residues on nodes, either outscores the charge by `thr` or lies
/// within `near` nt of the row, and scores what its place asks (bits_needed).
pub fn keep_piece(p: &Piece, o: &OrfOpts) -> bool {
    !p.repeat && p.margin >= o.revthr && p.nmatch >= o.minm && (p.bits - p.charge >= o.thr || p.dist <= o.near)
        && p.bits >= bits_needed(p)
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
const START_NT: i64 = 30;
const START_COVER: f64 = 0.7;
const STOP_NODES: i64 = 20;
const STOP_NT: i64 = 90;
const STOP_REACH_NODES: i64 = 100;
/// nt past a row's last exon read for a stop codon
pub const STOP_SCAN_NT: i64 = 1500;

/// nt upstream of a row's first exon to look for its start codon (the farthest
/// in-frame ATG before a stop), when the row earns the search: the exon starts
/// within START_NODES of the model's first node and the row's exons cover
/// START_COVER of the model. A lone hit or a thin row is too often no gene's
/// first exon. Alignments fade a few codons short of the start.
pub fn start_reach(e: &End) -> Option<i64> {
    (e.short < START_NODES && e.covered as f64 >= START_COVER * e.m as f64).then_some(START_NT)
}

/// Whether the first in-frame stop, `dist` nt past a row's last exon, ends the
/// row. Within STOP_NODES of the model's last node a stop within STOP_NT does.
/// Further from the model's end, up to STOP_REACH_NODES, any stop read does
/// unless a donor lies on the way (`donor`): a further exon. Beyond that the
/// exon is rarely the gene's last and is left as it is.
pub fn take_stop(e: &End, dist: i64, donor: impl FnOnce() -> bool) -> bool {
    (e.short <= STOP_NODES && dist <= STOP_NT) || (e.short <= STOP_REACH_NODES && !donor())
}
