//! Pseudogene disablements per gene copy. Every junction's annotate-mode path
//! reports frameshifts, in-frame stops and non-canonical introns. A disablement
//! is strict when no intron fits where it sits (inside an anchor exon, or in a
//! gap too short for an intron); a non-canonical intron always counts. A copy
//! is called, as in SAPPHYRE, on a strict frameshift or stop, or on at least
//! `min_strict` strict disablements.

use super::chain::Chain;
use super::splice::{SpliceResult, DIS_FS, DIS_STOP};

#[derive(Clone, Copy, Default, Debug)]
pub struct Pseudo {
    pub njx: i64,
    pub nfs: i64,
    pub fs_tight: i64,
    pub nstop: i64,
    pub stop_tight: i64,
    pub nnonc: i64,
}

impl Pseudo {
    pub fn count(res: &SpliceResult) -> Pseudo {
        let mut p = Pseudo { njx: 1, ..Default::default() };
        for d in &res.dis {
            match d.kind {
                DIS_FS => { p.nfs += 1; p.fs_tight += d.tight as i64; }
                DIS_STOP => { p.nstop += 1; p.stop_tight += d.tight as i64; }
                _ => p.nnonc += 1,
            }
        }
        p
    }
    pub fn merge(&mut self, o: &Pseudo) {
        self.nfs += o.nfs; self.fs_tight += o.fs_tight;
        self.nstop += o.nstop; self.stop_tight += o.stop_tight;
        self.nnonc += o.nnonc; self.njx += o.njx;
    }
    pub fn strict(&self) -> i64 { self.fs_tight + self.stop_tight + self.nnonc }
    pub fn called(&self, min_strict: i64) -> bool {
        self.fs_tight > 0 || self.stop_tight > 0 || self.strict() >= min_strict
    }
    /// gene scaffold strand n_junctions n_fs fs_strict n_stop stop_strict n_noncanon strict called flags
    pub fn row(&self, c: &Chain, min_strict: i64) -> String {
        let mut flags: Vec<String> = Vec::new();
        if self.fs_tight > 0 { flags.push(format!("fs:{}", self.fs_tight)); }
        if self.stop_tight > 0 { flags.push(format!("stop:{}", self.stop_tight)); }
        if self.nnonc > 0 { flags.push(format!("noncanonical:{}", self.nnonc)); }
        format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n", c.gene, c.scaffold, c.strand as char, self.njx,
                self.nfs, self.fs_tight, self.nstop, self.stop_tight, self.nnonc, self.strict(),
                self.called(min_strict) as i32, if flags.is_empty() { "-".into() } else { flags.join(";") })
    }
}
