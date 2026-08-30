//! Quality window cutting, and the pass/fail filters.
//!
//! Ported from fastp 1.3.6 `src/filter.cpp` (MIT, (c) 2016 OpenGene).
//!
//! Everything here needs a quality string, so all of it no-ops on FASTA input.

use crate::reads::read::ReadRec;
use crate::reads::seqops::{count_adjacent_diffs, count_quality_metrics};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FilterVerdict {
    Pass,
    FailLength,
    FailTooLong,
    FailQuality,
    FailNBase,
    FailComplexity,
    FailAdapterDimer,
}

impl FilterVerdict {
    /// fastp's severity scale from `common.h` — "bigger number means worse".
    /// A pair is recorded under the worse of its two mates' verdicts, so the
    /// exact numbers matter for reproducing its counts.
    pub fn severity(self) -> u8 {
        match self {
            FilterVerdict::Pass => 0,
            FilterVerdict::FailNBase => 12,
            FilterVerdict::FailLength => 16,
            FilterVerdict::FailTooLong => 17,
            FilterVerdict::FailQuality => 20,
            FilterVerdict::FailComplexity => 24,
            FilterVerdict::FailAdapterDimer => 28,
        }
    }

    /// The worse of two verdicts, fastp's `max(result1, result2)`.
    pub fn worse(self, other: Self) -> Self {
        if other.severity() > self.severity() {
            other
        } else {
            self
        }
    }
}

/// fastp's quality-filter block. Phred values, not raw bytes.
#[derive(Debug, Clone, Copy)]
pub struct QualFilter {
    pub enabled: bool,
    /// A base at or above this phred score counts as qualified.
    pub qualified_qual: u8,
    /// Percent of unqualified bases that fails the read.
    pub unqualified_percent_limit: f64,
    pub n_base_limit: usize,
    /// Mean phred below which the read fails. 0 disables the check.
    pub avg_qual_req: i64,
}

impl Default for QualFilter {
    fn default() -> Self {
        // fastp defaults: -q 15, -u 40, -n 5, -e 0
        QualFilter {
            enabled: true,
            qualified_qual: 15,
            unqualified_percent_limit: 40.0,
            n_base_limit: 5,
            avg_qual_req: 0,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct LengthFilter {
    pub enabled: bool,
    pub required_length: usize,
    /// 0 means no ceiling.
    pub max_length: usize,
}

impl Default for LengthFilter {
    fn default() -> Self {
        // fastp default: -l 15
        LengthFilter { enabled: true, required_length: 15, max_length: 0 }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct ComplexityFilter {
    pub enabled: bool,
    /// Minimum fraction of positions differing from the next base.
    pub threshold: f64,
}

impl Default for ComplexityFilter {
    fn default() -> Self {
        // fastp default: off, 0.3 when switched on
        ComplexityFilter { enabled: false, threshold: 0.3 }
    }
}

/// Sliding-window quality trimming, all three fastp modes.
#[derive(Debug, Clone, Copy, Default)]
pub struct QualityCut {
    pub enabled_front: bool,
    pub enabled_tail: bool,
    pub enabled_right: bool,
    pub window_size_front: usize,
    pub window_size_tail: usize,
    pub window_size_right: usize,
    pub quality_front: u8,
    pub quality_tail: u8,
    pub quality_right: u8,
}

#[derive(Debug, Clone, Copy, Default)]
pub struct FilterOptions {
    pub qual: Option<QualFilter>,
    pub length: Option<LengthFilter>,
    pub complexity: Option<ComplexityFilter>,
}

/// fastp's `passLowComplexityFilter`.
pub fn passes_complexity(r: &ReadRec<'_>, threshold: f64) -> bool {
    let len = r.len();
    if len <= 1 {
        return false;
    }
    let diff = count_adjacent_diffs(r.bases());
    (diff as f64) / ((len - 1) as f64) >= threshold
}

/// fastp's `passFilter`.
pub fn pass_filter(r: &ReadRec<'_>, opts: &FilterOptions) -> FilterVerdict {
    if r.is_empty() {
        return FilterVerdict::FailLength;
    }
    let rlen = r.len();

    if let Some(q) = opts.qual {
        if q.enabled && r.has_qual() {
            let m = count_quality_metrics(r.quals(), r.bases(), i64::from(q.qualified_qual) + 33);
            if (m.low_qual as f64) > q.unqualified_percent_limit * rlen as f64 / 100.0 {
                return FilterVerdict::FailQuality;
            }
            if q.avg_qual_req > 0 && m.total_qual / (rlen as i64) < q.avg_qual_req {
                return FilterVerdict::FailQuality;
            }
            if m.n_bases > q.n_base_limit {
                return FilterVerdict::FailNBase;
            }
        }
    }

    if let Some(l) = opts.length {
        if l.enabled {
            if rlen < l.required_length {
                return FilterVerdict::FailLength;
            }
            if l.max_length > 0 && rlen > l.max_length {
                return FilterVerdict::FailTooLong;
            }
        }
    }

    if let Some(c) = opts.complexity {
        if c.enabled && !passes_complexity(r, c.threshold) {
            return FilterVerdict::FailComplexity;
        }
    }

    FilterVerdict::Pass
}

/// Fixed-length trims plus quality window cutting, fastp's `trimAndCut`.
///
/// Returns false when the read should be dropped (fastp returns NULL).
pub fn trim_and_cut(
    r: &mut ReadRec<'_>,
    front: usize,
    tail: usize,
    cut: &QualityCut,
) -> bool {
    if front == 0 && tail == 0 && !cut.enabled_front && !cut.enabled_tail && !cut.enabled_right {
        return true;
    }

    let l = r.len();
    if front + tail > l {
        return false;
    }
    let mut rlen = l - front - tail;

    if !cut.enabled_front && !cut.enabled_tail && !cut.enabled_right {
        if front > 0 {
            r.trim_front(front);
        }
        r.resize(rlen);
        return true;
    }
    if !r.has_qual() {
        // No quality to cut on; the fixed trims still apply.
        if front > 0 {
            r.trim_front(front);
        }
        r.resize(rlen);
        return true;
    }

    let qual = r.quals();
    let seq = r.bases();
    let mut front = front;

    // Forward: advance until a window is good enough.
    if cut.enabled_front {
        let w = cut.window_size_front;
        if l < front + tail + w || w == 0 {
            return false;
        }
        let mut total: i64 = 0;
        for i in 0..w - 1 {
            total += i64::from(qual[front + i]);
        }
        let mut s = front;
        while s + w < l - tail {
            total += i64::from(qual[s + w - 1]);
            if s > front {
                total -= i64::from(qual[s - 1]);
            }
            if total >= (w as i64) * (33 + i64::from(cut.quality_front)) {
                break;
            }
            s += 1;
        }
        if s > 0 {
            s = s + w - 1;
        }
        while s < l && seq[s] == b'N' {
            s += 1;
        }
        front = s;
        rlen = l.saturating_sub(front + tail);
    }

    // Right: cut at the first bad window and keep the good bases inside it.
    if cut.enabled_right {
        let w = cut.window_size_right;
        if l < front + tail + w || w == 0 {
            return false;
        }
        let mut total: i64 = 0;
        for i in 0..w - 1 {
            total += i64::from(qual[front + i]);
        }
        let mut s = front;
        let mut found_low_qual_window = false;
        while s + w < l - tail {
            total += i64::from(qual[s + w - 1]);
            if s > front {
                total -= i64::from(qual[s - 1]);
            }
            if total < (w as i64) * (33 + i64::from(cut.quality_right)) {
                found_low_qual_window = true;
                break;
            }
            s += 1;
        }
        if found_low_qual_window {
            while s < l - 1 && i64::from(qual[s]) >= 33 + i64::from(cut.quality_right) {
                s += 1;
            }
            rlen = s - front;
        }
    }

    // Backward, only when right mode is off.
    if !cut.enabled_right && cut.enabled_tail {
        let w = cut.window_size_tail;
        if l < front + tail + w || w == 0 {
            return false;
        }
        let mut total: i64 = 0;
        let t0 = l - tail - 1;
        for i in 0..w - 1 {
            total += i64::from(qual[t0 - i]);
        }
        let mut t = t0 as isize;
        while t - (w as isize) >= front as isize {
            total += i64::from(qual[(t - w as isize + 1) as usize]);
            if (t as usize) < t0 {
                total -= i64::from(qual[(t + 1) as usize]);
            }
            if total >= (w as i64) * (33 + i64::from(cut.quality_tail)) {
                break;
            }
            t -= 1;
        }
        if (t as usize) < l - 1 {
            t = t - w as isize + 1;
        }
        while t >= 0 && seq[t as usize] == b'N' {
            t -= 1;
        }
        rlen = (t - front as isize + 1).max(0) as usize;
    }

    if rlen == 0 || front + 1 >= l {
        return false;
    }

    r.trim_front(front);
    r.resize(rlen);
    true
}

#[cfg(test)]
#[path = "../tests/filter.rs"]
mod tests;
