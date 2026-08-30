//! The stage order, and the options that drive it.
//!
//! Mirrors fastp 1.3.6 `src/seprocessor.cpp` / `src/peprocessor.cpp` (MIT,
//! (c) 2016 OpenGene) — the sequence of calls, not their threading, which is
//! prepare's business.
//!
//! Order per read: fixed trim and quality windows, polyG, adapter, polyX, then
//! the filters.
//!
//! A pair runs the same stages, with the adapter step replaced by overlap
//! analysis — sequence matching is only the fallback for pairs the overlap
//! could not place, which is fastp's `if(!trimmed)` branch. On a short insert
//! the overlap locates the adapter better than its sequence does, and needs no
//! adapter to be known at all. The dimer check follows, and fires only when
//! adapter was actually found in one of the mates.

use crate::reads::adapter::{trim_by_multi_sequences, trim_by_sequence, AdapterScratch};
use crate::reads::seed::PreparedAdapter;
use crate::reads::correct::{correct_by_overlap, PairBuffers};
use crate::reads::filter::{pass_filter, trim_and_cut, FilterOptions, FilterVerdict, QualityCut};
use crate::reads::overlap::{analyze, trim_by_overlap, OverlapResult, OverlapScratch};
use crate::reads::polyx::{trim_poly_g, trim_poly_x};
use crate::reads::read::ReadRec;

/// Everything the trimming chain can be told to do.
///
/// The default is fastp's default for a read library: adapters and polyG on,
/// quality and length filters on, polyX and complexity off.
#[derive(Debug, Clone)]
pub struct TrimOptions {
    /// Adapter for R1, and for unpaired reads. fastp's `adapter.sequence`.
    /// Seed tables are built once here, not per read.
    pub adapter_r1: Option<PreparedAdapter>,
    /// Adapter for R2. fastp's `adapter.sequenceR2`.
    pub adapter_r2: Option<PreparedAdapter>,
    /// Extra adapters applied to both mates on top of the above, fastp's
    /// `--adapter_fasta` list.
    pub adapters: Vec<PreparedAdapter>,
    pub adapter_trimming: bool,
    pub trim_poly_g: bool,
    pub poly_g_min_len: usize,
    pub trim_poly_x: bool,
    pub poly_x_min_len: usize,
    pub trim_front: usize,
    pub trim_tail: usize,
    pub quality_cut: QualityCut,
    pub filters: FilterOptions,
    /// Overlap analysis, only reachable when reads arrive as pairs.
    pub overlap: bool,
    pub overlap_diff_limit: usize,
    pub overlap_require: usize,
    pub overlap_diff_percent_limit: f64,
    pub allow_gap_overlap: bool,
    /// Rewrite disagreeing bases inside the overlap. Read only by
    /// `correct_pair`, which is part of the unwired merge feature.
    #[allow(dead_code)]
    pub correction: bool,
    /// Drop a pair whose mates are both this short after trimming. fastp's
    /// adapter-dimer check, new in 1.3.6.
    pub dimer_max_len: usize,
}

impl Default for TrimOptions {
    fn default() -> Self {
        TrimOptions {
            adapter_r1: None,
            adapter_r2: None,
            adapters: Vec::new(),
            adapter_trimming: true,
            trim_poly_g: true,
            poly_g_min_len: 10,
            trim_poly_x: false,
            poly_x_min_len: 10,
            trim_front: 0,
            trim_tail: 0,
            quality_cut: QualityCut::default(),
            filters: FilterOptions {
                qual: Some(Default::default()),
                length: Some(Default::default()),
                complexity: Some(Default::default()),
            },
            overlap: true,
            overlap_diff_limit: 5,
            overlap_require: 30,
            overlap_diff_percent_limit: 0.2,
            allow_gap_overlap: false,
            correction: false,
            dimer_max_len: 2,
        }
    }
}

/// Largest insert size tracked; anything longer, or any pair whose overlap
/// could not be placed, lands in the final bucket. fastp's `insertSizeMax`.
pub const INSERT_SIZE_MAX: usize = 512;

/// Tallies for the run summary.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct TrimStats {
    pub reads_in: u64,
    pub reads_passed: u64,
    pub adapter_trimmed_reads: u64,
    pub adapter_trimmed_bases: u64,
    pub polyx_trimmed_reads: u64,
    pub polyx_trimmed_bases: u64,
    pub failed_quality: u64,
    pub failed_n_base: u64,
    pub failed_length: u64,
    pub failed_too_long: u64,
    pub failed_complexity: u64,
    pub failed_adapter_dimer: u64,
    pub corrected_bases: u64,
    pub pairs_merged: u64,
    /// Bases seen before trimming, and bases kept after it.
    pub bases_in: u64,
    pub bases_out: u64,
    /// Length of every surviving read, indexed by length.
    pub length_hist: Vec<u64>,
    /// Fragment length per pair, indexed by size; the last bucket counts pairs
    /// whose overlap could not be placed. Only filled for paired input, and
    /// only when the overlap was computed at all.
    pub insert_hist: Vec<u64>,
}

impl TrimStats {
    fn note_kept(&mut self, len: usize) {
        self.bases_out += len as u64;
        if self.length_hist.len() <= len {
            self.length_hist.resize(len + 1, 0);
        }
        self.length_hist[len] += 1;
    }

    /// fastp's `statInsertSize`: a positive offset means the fragment spans
    /// both reads, a negative one means it is shorter than a single read.
    fn note_insert(&mut self, ov: &OverlapResult, len1: usize, len2: usize, f1: usize, f2: usize) {
        if self.insert_hist.is_empty() {
            self.insert_hist = vec![0; INSERT_SIZE_MAX + 1];
        }
        let isize = if !ov.overlapped {
            INSERT_SIZE_MAX
        } else if ov.offset > 0 {
            (len1 + len2 + f1 + f2).saturating_sub(ov.overlap_len)
        } else {
            ov.overlap_len + f1 + f2
        };
        self.insert_hist[isize.min(INSERT_SIZE_MAX)] += 1;
    }

    fn record(&mut self, verdict: FilterVerdict) -> bool {
        match verdict {
            FilterVerdict::Pass => {
                self.reads_passed += 1;
                true
            }
            FilterVerdict::FailQuality => {
                self.failed_quality += 1;
                false
            }
            FilterVerdict::FailNBase => {
                self.failed_n_base += 1;
                false
            }
            FilterVerdict::FailLength => {
                self.failed_length += 1;
                false
            }
            FilterVerdict::FailTooLong => {
                self.failed_too_long += 1;
                false
            }
            FilterVerdict::FailComplexity => {
                self.failed_complexity += 1;
                false
            }
            FilterVerdict::FailAdapterDimer => {
                self.failed_adapter_dimer += 1;
                false
            }
        }
    }

    /// Record one verdict `n` times, for fastp's pair accounting.
    fn record_n(&mut self, verdict: FilterVerdict, n: u64) -> bool {
        let mut pass = false;
        for _ in 0..n {
            pass = self.record(verdict);
        }
        pass
    }
}

/// Per-call scratch, so a long run allocates nothing per read.
#[derive(Default)]
pub struct TrimScratch {
    overlap: OverlapScratch,
    adapter: AdapterScratch,
}

/// The trim stages shared by single and paired reads, in fastp's order.
/// Returns false when the read was dropped before filtering.
fn trim_common(
    r: &mut ReadRec<'_>,
    opts: &TrimOptions,
    stats: &mut TrimStats,
    scratch: &mut TrimScratch,
) -> bool {
    if !trim_and_cut(r, opts.trim_front, opts.trim_tail, &opts.quality_cut) {
        stats.failed_length += 1;
        return false;
    }
    if opts.trim_poly_g {
        let removed = trim_poly_g(r, opts.poly_g_min_len);
        if removed > 0 {
            stats.polyx_trimmed_reads += 1;
            stats.polyx_trimmed_bases += removed as u64;
        }
    }
    if opts.adapter_trimming {
        let mut bases = 0u64;
        if let Some(a) = &opts.adapter_r1 {
            if let Some(h) = trim_by_sequence(r, a, 4, &mut scratch.adapter) {
                bases += h.removed as u64;
            }
        }
        bases += trim_by_fasta_list(r, opts, &mut scratch.adapter);
        if bases > 0 {
            stats.adapter_trimmed_reads += 1;
            stats.adapter_trimmed_bases += bases;
        }
    }
    if opts.trim_poly_x {
        if let Some((_, removed)) = trim_poly_x(r, opts.poly_x_min_len) {
            stats.polyx_trimmed_reads += 1;
            stats.polyx_trimmed_bases += removed as u64;
        }
    }
    true
}

/// The `--adapter_fasta` list, applied on top of the per-mate adapter.
fn trim_by_fasta_list(
    r: &mut ReadRec<'_>,
    opts: &TrimOptions,
    scratch: &mut AdapterScratch,
) -> u64 {
    if opts.adapters.is_empty() {
        return 0;
    }
    trim_by_multi_sequences(r, &opts.adapters, scratch) as u64
}

/// Process one unpaired read. `Some(window)` is what should be kept.
pub fn process_single<'a>(
    seq: &'a [u8],
    qual: &'a [u8],
    opts: &TrimOptions,
    stats: &mut TrimStats,
    scratch: &mut TrimScratch,
) -> Option<ReadRec<'a>> {
    stats.reads_in += 1;
    stats.bases_in += seq.len() as u64;
    let mut r = ReadRec::new(seq, qual);
    if !trim_common(&mut r, opts, stats, scratch) {
        return None;
    }
    if stats.record(pass_filter(&r, &opts.filters)) {
        stats.note_kept(r.len());
        Some(r)
    } else {
        None
    }
}

/// What survived a pair.
pub struct PairOutcome<'a> {
    pub r1: Option<ReadRec<'a>>,
    pub r2: Option<ReadRec<'a>>,
}

/// Process a mate pair, in fastp's `peprocessor` order.
pub fn process_pair<'a>(
    seq1: &'a [u8],
    qual1: &'a [u8],
    seq2: &'a [u8],
    qual2: &'a [u8],
    opts: &TrimOptions,
    stats: &mut TrimStats,
    scratch: &mut TrimScratch,
) -> PairOutcome<'a> {
    stats.reads_in += 2;
    stats.bases_in += (seq1.len() + seq2.len()) as u64;
    let mut r1 = ReadRec::new(seq1, qual1);
    let mut r2 = ReadRec::new(seq2, qual2);

    // 1. Fixed trims and quality windows.
    if !trim_and_cut(&mut r1, opts.trim_front, opts.trim_tail, &opts.quality_cut) {
        stats.failed_length += 1;
        return PairOutcome { r1: None, r2: None };
    }
    if !trim_and_cut(&mut r2, opts.trim_front, opts.trim_tail, &opts.quality_cut) {
        stats.failed_length += 1;
        return PairOutcome { r1: None, r2: None };
    }
    let front1 = r1.front_trimmed();
    let front2 = r2.front_trimmed();

    // 2. polyG, before the overlap so the analysis sees the real fragment.
    if opts.trim_poly_g {
        for r in [&mut r1, &mut r2] {
            let removed = trim_poly_g(r, opts.poly_g_min_len);
            if removed > 0 {
                stats.polyx_trimmed_reads += 1;
                stats.polyx_trimmed_bases += removed as u64;
            }
        }
    }

    // 3. Adapter: overlap first, sequence matching only as the fallback.
    let mut trimmed1 = false;
    let mut trimmed2 = false;
    if opts.adapter_trimming {
        let mut by_overlap = false;
        if opts.overlap {
            let ov = analyze(
                r1.bases(),
                r2.bases(),
                opts.overlap_diff_limit,
                opts.overlap_require,
                opts.overlap_diff_percent_limit,
                opts.allow_gap_overlap,
                &mut scratch.overlap,
            );
            stats.note_insert(&ov, r1.len(), r2.len(), front1, front2);
            if let Some((removed1, removed2)) =
                trim_by_overlap(&mut r1, &mut r2, ov, front1, front2)
            {
                by_overlap = true;
                trimmed1 = true;
                trimmed2 = true;
                stats.adapter_trimmed_reads += 2;
                stats.adapter_trimmed_bases += (removed1 + removed2) as u64;
            }
        }
        if !by_overlap {
            if let Some(a) = &opts.adapter_r1 {
                if let Some(h) = trim_by_sequence(&mut r1, a, 4, &mut scratch.adapter) {
                    trimmed1 = true;
                    stats.adapter_trimmed_reads += 1;
                    stats.adapter_trimmed_bases += h.removed as u64;
                }
            }
            if let Some(a) = &opts.adapter_r2 {
                if let Some(h) = trim_by_sequence(&mut r2, a, 4, &mut scratch.adapter) {
                    trimmed2 = true;
                    stats.adapter_trimmed_reads += 1;
                    stats.adapter_trimmed_bases += h.removed as u64;
                }
            }
        }
        let extra1 = trim_by_fasta_list(&mut r1, opts, &mut scratch.adapter);
        let extra2 = trim_by_fasta_list(&mut r2, opts, &mut scratch.adapter);
        if extra1 > 0 {
            trimmed1 = true;
            stats.adapter_trimmed_bases += extra1;
        }
        if extra2 > 0 {
            trimmed2 = true;
            stats.adapter_trimmed_bases += extra2;
        }
    }

    // 4. polyX.
    if opts.trim_poly_x {
        for r in [&mut r1, &mut r2] {
            if let Some((_, removed)) = trim_poly_x(r, opts.poly_x_min_len) {
                stats.polyx_trimmed_reads += 1;
                stats.polyx_trimmed_bases += removed as u64;
            }
        }
    }

    // 5. Adapter dimer: only when adapter was actually found, so a pair of
    //    genuinely tiny inserts is not mistaken for one.
    let is_dimer = (trimmed1 || trimmed2)
        && r1.len() <= opts.dimer_max_len
        && r2.len() <= opts.dimer_max_len;

    // 6. Filters. fastp judges the pair, not the mate: it records the worse of
    //    the two verdicts against both reads, and keeps the pair only if both
    //    pass. Counting per mate instead inflates the pass count and halves
    //    every failure tally.
    let (v1, v2) = if is_dimer {
        (FilterVerdict::FailAdapterDimer, FilterVerdict::FailAdapterDimer)
    } else {
        (pass_filter(&r1, &opts.filters), pass_filter(&r2, &opts.filters))
    };
    let worse = v1.worse(v2);
    let keep = stats.record_n(worse, 2);

    if keep {
        stats.note_kept(r1.len());
        stats.note_kept(r2.len());
        PairOutcome { r1: Some(r1), r2: Some(r2) }
    } else {
        PairOutcome { r1: None, r2: None }
    }
}

/// Correct disagreeing bases in a pair, for callers that want it. Kept apart
/// from [`process_pair`] because it has to own its buffers. Part of the
/// unwired merge feature; see `correct.rs`.
#[allow(dead_code)]
pub fn correct_pair(p: &mut PairBuffers, opts: &TrimOptions, stats: &mut TrimStats) {
    if !opts.correction {
        return;
    }
    let mut scratch = OverlapScratch::default();
    let ov = analyze(
        &p.seq1,
        &p.seq2,
        opts.overlap_diff_limit,
        opts.overlap_require,
        opts.overlap_diff_percent_limit,
        opts.allow_gap_overlap,
        &mut scratch,
    );
    let res = correct_by_overlap(p, ov);
    stats.corrected_bases += res.corrected as u64;
}

#[cfg(test)]
#[path = "../tests/pipeline.rs"]
mod tests;
