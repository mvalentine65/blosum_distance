//! One junction between anchor exons A and B (A before B in model order):
//! orient to the gene strand, pick each anchor's frame, trim the anchors to
//! the `keep` nodes next to the junction, and run the spliced alignment.

use super::hmm::Hmm;
use super::orf::{anchor_frame, Aligner};
use super::sites::revcomp;
use super::splice::{splice, Locus, Params, SpliceResult, SpliceStatus};

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum JxStatus { Ok, Nodes, Order, Short, Big, NoPath }

impl JxStatus {
    pub fn name(self) -> &'static str {
        match self {
            JxStatus::Ok => "ok", JxStatus::Nodes => "nodes", JxStatus::Order => "order",
            JxStatus::Short => "short", JxStatus::Big => "big", JxStatus::NoPath => "nopath",
        }
    }
}

pub struct Junction {
    pub status: JxStatus,
    pub s: Vec<u8>, // window, oriented to the gene strand
    pub w: i64,
    pub lo: i64, // locus start in s
    pub strand: u8,
    pub gap: i64,
    pub d: i64, // locus length
    pub axe: i64, // A's anchor end and B's anchor start, locus coords
    pub bxs: i64,
    pub res: SpliceResult,
}

impl Junction {
    /// Window + strand coordinate (1-based) of locus position x (0-based).
    pub fn pos(&self, x: i64) -> i64 {
        let o = self.lo + x + 1;
        if self.strand == b'-' { self.w - o + 1 } else { o }
    }
    pub fn dna(&self) -> &[u8] {
        &self.s[self.lo as usize..(self.lo + self.d) as usize]
    }
}

/// seq: + strand window; A and B nodes and window + strand extents (1-based).
#[allow(clippy::too_many_arguments)]
pub fn junction(al: &mut Aligner, hmm_id: usize, hmm: &Hmm, prm: &Params, keep: usize, seq: &[u8], strand: u8,
                ak1: i64, mut ak2: i64, mut alo: i64, mut ahi: i64, mut bk1: i64, bk2: i64, mut blo: i64, mut bhi: i64) -> Junction {
    let w = seq.len() as i64;
    let mut jx = Junction { status: JxStatus::Ok, s: Vec::new(), w, lo: 0, strand, gap: 0, d: 0, axe: 0, bxs: 0, res: SpliceResult::default() };
    // hits either side of a frameshift often share a node or two: split the overlap
    if ak2 >= bk1 { let m = (ak2 + bk1) / 2; ak2 = m; bk1 = m + 1; }
    if ak1 < 1 || bk2 > hmm.m as i64 || ak1 > ak2 || bk1 > bk2 { jx.status = JxStatus::Nodes; return jx; }
    jx.s = if strand == b'-' {
        let t = w - ahi + 1; ahi = w - alo + 1; alo = t;
        let t = w - bhi + 1; bhi = w - blo + 1; blo = t;
        revcomp(seq)
    } else {
        seq.to_vec()
    };
    // ... and a few bases of genome
    if ahi >= blo && alo < blo && ahi < bhi { let m = (ahi + blo) / 2; ahi = m; blo = m + 1; }
    if ahi >= blo || alo > ahi || blo > bhi { jx.status = JxStatus::Order; return jx; }

    let (afo, a_lo, _a_hi, a_k) = anchor_frame(al, hmm_id, hmm, &jx.s[(alo - 1) as usize..ahi as usize], ak1 as usize, ak2 as usize, 0, keep, 30);
    let (bfo, _b_lo, b_hi, b_k) = anchor_frame(al, hmm_id, hmm, &jx.s[(blo - 1) as usize..bhi as usize], bk1 as usize, bk2 as usize, 1, keep, 30);

    jx.lo = alo - 1 + a_lo;
    jx.d = (blo - 1 + b_hi) - jx.lo;
    let axe = ahi - jx.lo;
    let bxs = blo - 1 - jx.lo;
    jx.gap = bxs - axe;
    jx.axe = axe;
    jx.bxs = bxs;
    if jx.gap < prm.min_gap_nt { jx.status = JxStatus::Short; return jx; }
    let loc = Locus {
        dna: &jx.s[jx.lo as usize..(jx.lo + jx.d) as usize],
        axs: 0, axe, bxs, bxe: jx.d, afo, bfo,
        k1: a_k, k2: b_k, ak2: ak2 as usize, bk1: bk1 as usize,
    };
    let mut res = SpliceResult::default();
    jx.status = match splice(hmm, &loc, prm, &mut res, &mut al.sw) {
        SpliceStatus::Ok => JxStatus::Ok,
        SpliceStatus::Big => JxStatus::Big,
        SpliceStatus::NoPath => JxStatus::NoPath,
    };
    jx.res = res;
    jx
}
