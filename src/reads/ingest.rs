//! Driving the trim chain over a whole input, in one pass.
//!
//! This is the piece that makes trimming nearly free: `dedupe` already streams
//! every read once, so the chain rides along on that pass rather than adding
//! its own. The only wrinkle is ordering — adapters have to be *detected*
//! before the first read can be trimmed, and a read has to be trimmed before it
//! is hashed, because trimming is what makes adapter-bearing duplicates
//! collapse onto each other.
//!
//! fastp resolves that by opening the file twice: an evaluator pass reads the
//! head to find adapters, then the real pass re-reads from the top. We own the
//! loop, so instead we **stage and backfill**: hold the first
//! [`DETECT_SAMPLE`] records without hashing them, detect from that buffer,
//! then drain the buffer through the chain and stream the remainder inline.
//! One decompression, no reopen, no seek.
//!
//! Nothing here touches the FASTQ readers or pyo3 — a caller supplies records
//! and consumes trimmed bases through a sink, which keeps the whole thing
//! testable without the crate's heavier dependencies.

use std::path::PathBuf;

use crate::reads::detect::{detect_adapter, Detected};
use crate::reads::pipeline::{
    process_pair, process_single, TrimOptions, TrimScratch, TrimStats,
};
use crate::reads::seed::PreparedAdapter;

/// Records held back for detection before any hashing happens. Matches fastp's
/// evaluator sample size; at 150 bp both mates cost roughly 80 MB, against a
/// dedup arena already sized in the hundreds of megabytes.
pub const DETECT_SAMPLE: usize = 256 * 1024;

/// How one taxa's files pair up.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum InputGroup {
    /// R1 and R2 of the same library, to be read in lockstep.
    Paired(PathBuf, PathBuf),
    /// A lone file: single-end reads, or a FASTA.
    Single(PathBuf),
}

/// Strip the mate marker from a file name, returning `(stem, mate, segment)`.
///
/// Mirrors the convention `prepare.py`'s `_TRUNCATE_TAXA_RE` already relies on
/// to group both mates under one taxa: a trailing `_1`/`_2` or `_R1`/`_R2`,
/// optionally followed by the segment index bcl2fastq appends to every file
/// ("_001", bumped per chunk once a library outgrows one file). The segment
/// comes back so pairing can require it to match -- `R1_002` against `R2_003`
/// would feed the analyser two unrelated fragments in lockstep.
fn split_mate(name: &str) -> Option<(String, u8, Option<String>)> {
    let mut stem = name;
    // Peel the extensions prepare accepts, longest suffix first.
    loop {
        let bytes = stem.as_bytes();
        let stripped = [
            ".fq.gz", ".fastq.gz", ".fa.gz", ".fasta.gz", ".fna.gz", ".fas.gz",
            ".fq", ".fastq", ".fa", ".fasta", ".fna", ".fas",
        ]
        .iter()
        .find(|ext| {
            bytes.len() >= ext.len()
                && bytes[bytes.len() - ext.len()..].eq_ignore_ascii_case(ext.as_bytes())
        })
        .map(|ext| stem.len() - ext.len());
        match stripped {
            Some(cut) => stem = &stem[..cut],
            None => break,
        }
    }
    // A bare marker wins over a segmented one, so `lib_1` stays mate 1 of `lib`
    // rather than becoming segment 1 of a stem with no marker left.
    if let Some((base, mate)) = trailing_mate(stem) {
        return Some((base, mate, None));
    }
    let (head, seg) = trailing_segment(stem)?;
    let (base, mate) = trailing_mate(head)?;
    Some((base, mate, Some(seg)))
}

/// Split a trailing `_1`/`_2`/`_R1`/`_R2` off a stem.
fn trailing_mate(stem: &str) -> Option<(String, u8)> {
    let bytes = stem.as_bytes();
    if bytes.len() >= 2 && (bytes[bytes.len() - 1] == b'1' || bytes[bytes.len() - 1] == b'2') {
        let mate = bytes[bytes.len() - 1] - b'0';
        let head = &stem[..stem.len() - 1];
        if let Some(base) = head.strip_suffix('_') {
            return Some((base.to_string(), mate));
        }
        if let Some(base) = head.strip_suffix("_R").or_else(|| head.strip_suffix("_r")) {
            return Some((base.to_string(), mate));
        }
    }
    None
}

/// Split a trailing `_<digits>` segment index off a stem.
fn trailing_segment(stem: &str) -> Option<(&str, String)> {
    let cut = stem.rfind('_')?;
    let seg = &stem[cut + 1..];
    if seg.is_empty() || !seg.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    Some((&stem[..cut], seg.to_string()))
}

/// Pair up a taxa's files by name, leaving anything unmatched on its own.
///
/// Order is preserved so a run is reproducible: groups come back in the order
/// their first file appeared.
pub fn group_inputs(paths: &[PathBuf]) -> Vec<InputGroup> {
    let mut out: Vec<InputGroup> = Vec::new();
    let mut used = vec![false; paths.len()];

    for i in 0..paths.len() {
        if used[i] {
            continue;
        }
        let name_i = paths[i]
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or_default();
        let Some((stem_i, mate_i, seg_i)) = split_mate(name_i) else {
            used[i] = true;
            out.push(InputGroup::Single(paths[i].clone()));
            continue;
        };

        // Look for the opposite mate of the same stem and segment in the same
        // directory.
        let mut partner = None;
        for (j, item) in paths.iter().enumerate().skip(i + 1) {
            if used[j] {
                continue;
            }
            if item.parent() != paths[i].parent() {
                continue;
            }
            let name_j = item.file_name().and_then(|s| s.to_str()).unwrap_or_default();
            if let Some((stem_j, mate_j, seg_j)) = split_mate(name_j) {
                if stem_j == stem_i && seg_j == seg_i && mate_j != mate_i {
                    partner = Some(j);
                    break;
                }
            }
        }

        match partner {
            Some(j) => {
                used[i] = true;
                used[j] = true;
                let (a, b) = if mate_i == 1 {
                    (paths[i].clone(), paths[j].clone())
                } else {
                    (paths[j].clone(), paths[i].clone())
                };
                out.push(InputGroup::Paired(a, b));
            }
            None => {
                used[i] = true;
                out.push(InputGroup::Single(paths[i].clone()));
            }
        }
    }
    out
}

struct StagedPair {
    s1: Vec<u8>,
    q1: Vec<u8>,
    s2: Vec<u8>,
    q2: Vec<u8>,
}

/// Runs the trim chain over a stream of records, staging the head so adapters
/// can be detected before anything is hashed.
pub struct Ingest {
    opts: TrimOptions,
    stats: TrimStats,
    scratch: TrimScratch,
    staged: Vec<StagedPair>,
    /// Set once detection has run and the staged backlog has been drained.
    armed: bool,
    detected_r1: Option<Vec<u8>>,
    detected_r2: Option<Vec<u8>>,
}

impl Ingest {
    pub fn new(opts: TrimOptions) -> Self {
        Ingest {
            opts,
            stats: TrimStats::default(),
            scratch: TrimScratch::default(),
            staged: Vec::new(),
            armed: false,
            detected_r1: None,
            detected_r2: None,
        }
    }

    pub fn stats(&self) -> &TrimStats {
        &self.stats
    }

    /// Adapters detection settled on, for the run summary.
    pub fn detected(&self) -> (Option<&[u8]>, Option<&[u8]>) {
        (self.detected_r1.as_deref(), self.detected_r2.as_deref())
    }

    /// Feed one mate pair. Survivors reach `sink` in input order.
    pub fn push_pair(
        &mut self,
        s1: &[u8],
        q1: &[u8],
        s2: &[u8],
        q2: &[u8],
        sink: &mut impl FnMut(&[u8]),
    ) {
        if self.armed {
            self.run_pair(s1, q1, s2, q2, sink);
            return;
        }
        self.staged.push(StagedPair {
            s1: s1.to_vec(),
            q1: q1.to_vec(),
            s2: s2.to_vec(),
            q2: q2.to_vec(),
        });
        if self.staged.len() >= DETECT_SAMPLE {
            self.arm(sink);
        }
    }

    /// Feed one unpaired read. Staged and detected the same way.
    pub fn push_single(&mut self, s: &[u8], q: &[u8], sink: &mut impl FnMut(&[u8])) {
        if self.armed {
            self.run_single(s, q, sink);
            return;
        }
        self.staged.push(StagedPair {
            s1: s.to_vec(),
            q1: q.to_vec(),
            s2: Vec::new(),
            q2: Vec::new(),
        });
        if self.staged.len() >= DETECT_SAMPLE {
            self.arm(sink);
        }
    }

    /// Flush anything still staged. Call once the input is exhausted.
    pub fn finish(&mut self, sink: &mut impl FnMut(&[u8])) {
        if !self.armed {
            self.arm(sink);
        }
    }

    /// Detect from the staged head, then drain it through the chain.
    fn arm(&mut self, sink: &mut impl FnMut(&[u8])) {
        let staged = std::mem::take(&mut self.staged);
        let paired = staged.first().is_some_and(|p| !p.s2.is_empty());

        if self.opts.adapter_trimming {
            let sample1: Vec<&[u8]> = staged.iter().map(|p| p.s1.as_slice()).collect();
            if let Detected::Known { seq, .. } | Detected::Novel { seq } =
                detect_adapter(&sample1, 1)
            {
                self.opts.adapter_r1 = Some(PreparedAdapter::new(&seq));
                self.detected_r1 = Some(seq);
            }
            if paired {
                let sample2: Vec<&[u8]> = staged.iter().map(|p| p.s2.as_slice()).collect();
                if let Detected::Known { seq, .. } | Detected::Novel { seq } =
                    detect_adapter(&sample2, 1)
                {
                    self.opts.adapter_r2 = Some(PreparedAdapter::new(&seq));
                    self.detected_r2 = Some(seq);
                }
            }
        }
        self.armed = true;

        for p in &staged {
            if p.s2.is_empty() {
                self.run_single(&p.s1, &p.q1, sink);
            } else {
                self.run_pair(&p.s1, &p.q1, &p.s2, &p.q2, sink);
            }
        }
    }

    fn run_pair(
        &mut self,
        s1: &[u8],
        q1: &[u8],
        s2: &[u8],
        q2: &[u8],
        sink: &mut impl FnMut(&[u8]),
    ) {
        let out = process_pair(
            s1,
            q1,
            s2,
            q2,
            &self.opts,
            &mut self.stats,
            &mut self.scratch,
        );
        if let Some(r) = out.r1 {
            sink(r.bases());
        }
        if let Some(r) = out.r2 {
            sink(r.bases());
        }
    }

    fn run_single(&mut self, s: &[u8], q: &[u8], sink: &mut impl FnMut(&[u8])) {
        if let Some(r) = process_single(s, q, &self.opts, &mut self.stats, &mut self.scratch) {
            sink(r.bases());
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reads::seqops::reverse_complement;

    fn p(s: &str) -> PathBuf {
        PathBuf::from(s)
    }

    #[test]
    fn pairs_the_real_library_naming() {
        let paths = vec![
            p("/in/DDM8637_CKDL230042833-1A_HGC33DSX7_L2_1.fq.gz"),
            p("/in/DDM8637_CKDL230042833-1A_HGC33DSX7_L2_2.fq.gz"),
        ];
        assert_eq!(
            group_inputs(&paths),
            vec![InputGroup::Paired(paths[0].clone(), paths[1].clone())]
        );
    }

    #[test]
    fn pairs_r1_r2_and_orders_them() {
        // R2 listed first: the group must still come back R1-then-R2.
        let paths = vec![p("/in/lib_R2.fastq"), p("/in/lib_R1.fastq")];
        assert_eq!(
            group_inputs(&paths),
            vec![InputGroup::Paired(p("/in/lib_R1.fastq"), p("/in/lib_R2.fastq"))]
        );
    }

    #[test]
    fn an_unmatched_file_stays_single() {
        let paths = vec![p("/in/lib_1.fq"), p("/in/other.fa")];
        assert_eq!(
            group_inputs(&paths),
            vec![
                InputGroup::Single(p("/in/lib_1.fq")),
                InputGroup::Single(p("/in/other.fa")),
            ]
        );
    }

    #[test]
    fn same_stem_in_different_directories_does_not_pair() {
        let paths = vec![p("/a/lib_1.fq"), p("/b/lib_2.fq")];
        assert!(matches!(group_inputs(&paths)[0], InputGroup::Single(_)));
        assert_eq!(group_inputs(&paths).len(), 2);
    }

    #[test]
    fn several_libraries_pair_independently() {
        let paths = vec![
            p("/in/a_1.fq.gz"),
            p("/in/b_1.fq.gz"),
            p("/in/a_2.fq.gz"),
            p("/in/b_2.fq.gz"),
        ];
        let groups = group_inputs(&paths);
        assert_eq!(groups.len(), 2);
        assert_eq!(
            groups[0],
            InputGroup::Paired(p("/in/a_1.fq.gz"), p("/in/a_2.fq.gz"))
        );
        assert_eq!(
            groups[1],
            InputGroup::Paired(p("/in/b_1.fq.gz"), p("/in/b_2.fq.gz"))
        );
    }

    #[test]
    fn a_digit_that_is_not_a_mate_marker_is_left_alone() {
        // No underscore before the digit, so this is not a mate suffix.
        assert_eq!(split_mate("sample2.fq"), None);
        assert_eq!(split_mate("lib_L001.fq"), None);
    }

    /// The staged head must reach the sink too — losing it would silently drop
    /// the first quarter-million reads of every library.
    #[test]
    fn short_input_still_drains_through_finish() {
        let insert = b"ACGTTGCAAGGTCCATTGACGATCGGATCCAAGGTTCCAAGGTTCCAAGGTTCCAA";
        let r2 = reverse_complement(insert);
        let q = vec![b'I'; insert.len()];

        let mut ing = Ingest::new(TrimOptions::default());
        let mut kept: Vec<Vec<u8>> = Vec::new();
        {
            let mut sink = |b: &[u8]| kept.push(b.to_vec());
            for _ in 0..10 {
                ing.push_pair(insert, &q, &r2, &q, &mut sink);
            }
            // Nothing has been emitted yet: it is all staged awaiting detection.
            ing.finish(&mut sink);
        }
        assert_eq!(kept.len(), 20, "10 pairs, both mates surviving");
        assert_eq!(ing.stats().reads_in, 20);
    }

    #[test]
    fn finish_is_idempotent() {
        let mut ing = Ingest::new(TrimOptions::default());
        let mut kept = 0usize;
        {
            let mut sink = |_: &[u8]| kept += 1;
            let q = vec![b'I'; 8];
            ing.push_single(b"ACGTACGT", &q, &mut sink);
            ing.finish(&mut sink);
            ing.finish(&mut sink);
        }
        // The read is below the length filter, so nothing is emitted -- the
        // point is that the second finish does not replay the first.
        assert_eq!(ing.stats().reads_in, 1);
        assert_eq!(kept, 0);
    }

    /// Adapter-bearing reads must collapse once trimmed, which is the whole
    /// reason trimming has to happen before hashing.
    #[test]
    fn trimming_makes_two_reads_identical() {
        let insert = b"ACGTTGCAAGGTCCATTGACGATCGGATCCAAGGTTCCAAGGTTCCAAGGTTCCAA";
        let mut a = insert.to_vec();
        a.extend_from_slice(b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA");
        let mut b = insert.to_vec();
        b.extend_from_slice(b"AGATCGGAAGAGCAC"); // same adapter, shorter remnant
        let qa = vec![b'I'; a.len()];
        let qb = vec![b'I'; b.len()];

        let opts = TrimOptions {
            adapter_r1: Some(PreparedAdapter::new(b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA")),
            ..Default::default()
        };
        let mut ing = Ingest::new(opts);
        let mut kept: Vec<Vec<u8>> = Vec::new();
        {
            let mut sink = |x: &[u8]| kept.push(x.to_vec());
            ing.push_single(&a, &qa, &mut sink);
            ing.push_single(&b, &qb, &mut sink);
            ing.finish(&mut sink);
        }
        assert_eq!(kept.len(), 2);
        assert_eq!(kept[0], kept[1], "untrimmed these would be distinct");
        assert_eq!(kept[0], insert.to_vec());
    }

}
