//! Adapter trimming by sequence.
//!
//! Ported from fastp 1.3.6 `src/adaptertrimmer.cpp` (MIT, (c) 2016 OpenGene).
//!
//! One ungapped pass, allowing one mismatch per eight compared bases. The scan
//! starts at a small negative offset because A-tailing commonly eats the
//! adapter's leading A.
//!
//! fastp follows this with two fallbacks that permit a single inserted or
//! deleted base. Those are **deliberately not ported**. They hand the matcher
//! the read and the adapter both from index 0, so they can only fire when a
//! read is essentially all adapter — which the dimer filter already catches —
//! and on this library they fired 288 times per 3 M reads, every one a false
//! positive: the removed pieces averaged Q35 and matched the adapter at 1-3
//! bases out of 8-10, i.e. no better than chance. Dropping them cut 26% of the
//! runtime and reproduced fastp 0.23.4's output exactly.


/// Scratch for the seed-filtered scan.
#[derive(Default)]
pub struct AdapterScratch {
    candidates: Vec<usize>,
    seen: Vec<u64>,
}
use crate::reads::read::ReadRec;
use crate::reads::seed::PreparedAdapter;
use crate::reads::seqops::count_mismatches_bounded;

const ALLOW_ONE_MISMATCH_FOR_EACH: usize = 8;

/// What a trim removed, for the run summary.
///
/// Only the count is ever reported, so this carries a length rather than the
/// bases: copying every trimmed tail cost one allocation per hit, 3.7 M of
/// them on the test library.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AdapterHit {
    pub removed: usize,
}

/// fastp scales the shortest acceptable match to the size of the adapter list:
/// more candidates mean more chances for a short spurious hit.
pub fn match_req_for(list_len: usize) -> usize {
    if list_len > 256 {
        6
    } else if list_len > 16 {
        5
    } else {
        4
    }
}

/// Trim one adapter off the 3' end of `r`.
///
/// Returns the removed bases when a match was found.
pub fn trim_by_sequence(
    r: &mut ReadRec<'_>,
    prepared: &PreparedAdapter,
    match_req: usize,
    scratch: &mut AdapterScratch,
) -> Option<AdapterHit> {
    let adapter = &prepared.seq[..];
    let rlen = r.len();
    let alen = adapter.len();
    if alen < match_req || rlen == 0 {
        return None;
    }

    let rdata = r.bases();

    // A-tailing usually swallows the adapter's first base, so fastp begins the
    // scan before the read starts and clips the comparison to what overlaps.
    let start: isize = if alen >= 16 {
        -4
    } else if alen >= 12 {
        -3
    } else if alen >= 8 {
        -2
    } else {
        0
    };

    let mut found_pos: Option<isize> = None;

    // Pass 1: hamming distance, no indels.
    //
    // Split in three. The negative-offset head and the tail where the adapter
    // hangs off the read both have a shrinking `cmplen`, and so a shrinking
    // mismatch budget the seed blocks were not cut for — those stay
    // exhaustive, and together they are only a few dozen offsets. The bulk in
    // between, where `cmplen == alen`, goes through the seed filter.
    let seeded_max: isize = rlen as isize - alen as isize;
    let scan_limit = rlen as isize - match_req as isize;

    let probe = |pos: isize| -> bool {
        let cmplen = ((rlen as isize - pos) as usize).min(alen);
        let allowed_mismatch = cmplen / ALLOW_ONE_MISMATCH_FOR_EACH;
        let start_offset = (-pos).max(0) as usize;
        if start_offset >= cmplen {
            return false;
        }
        let a = &adapter[start_offset..];
        let b = &rdata[(start_offset as isize + pos) as usize..];
        count_mismatches_bounded(a, b, cmplen - start_offset, allowed_mismatch) <= allowed_mismatch
    };

    // Head: negative offsets.
    let mut pos = start;
    while pos < 0.min(scan_limit) {
        if probe(pos) {
            found_pos = Some(pos);
            break;
        }
        pos += 1;
    }

    // Body: full-length comparisons, seed filtered when the adapter allows it.
    if found_pos.is_none() && seeded_max >= 0 {
        let body_max = seeded_max.min(scan_limit - 1);
        if body_max >= 0 {
            if prepared.can_filter() {
                prepared.candidates(rdata, body_max as usize, &mut scratch.seen, &mut scratch.candidates);
                let candidates = std::mem::take(&mut scratch.candidates);
                for &cand in &candidates {
                    if probe(cand as isize) {
                        found_pos = Some(cand as isize);
                        break;
                    }
                }
                scratch.candidates = candidates;
            } else {
                let mut pos = 0isize;
                while pos <= body_max {
                    if probe(pos) {
                        found_pos = Some(pos);
                        break;
                    }
                    pos += 1;
                }
            }
        }
    }

    // Tail: the adapter runs off the end, so the budget shrinks with it.
    if found_pos.is_none() {
        let mut pos = (seeded_max + 1).max(0);
        while pos < scan_limit {
            if probe(pos) {
                found_pos = Some(pos);
                break;
            }
            pos += 1;
        }
    }

    let pos = found_pos?;
    if pos < 0 {
        // The adapter starts before the read: the whole record is adapter.
        let removed = (alen as isize + pos) as usize;
        r.resize(0);
        Some(AdapterHit { removed })
    } else {
        let removed = rlen - pos as usize;
        r.resize(pos as usize);
        Some(AdapterHit { removed })
    }
}

/// Try every adapter in the list, fastp's `trimByMultiSequences`.
pub fn trim_by_multi_sequences(
    r: &mut ReadRec<'_>,
    adapters: &[PreparedAdapter],
    scratch: &mut AdapterScratch,
) -> usize {
    let match_req = match_req_for(adapters.len());
    let mut removed = 0;
    for adapter in adapters {
        if let Some(hit) = trim_by_sequence(r, adapter, match_req, scratch) {
            removed += hit.removed;
        }
    }
    removed
}

#[cfg(test)]
mod tests {
    use super::*;

    /// fastp's own `AdapterTrimmer::test()`, first case.
    #[test]
    fn fastp_trim_by_sequence_case() {
        let seq = b"TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAATTTTCCCCGGGG";
        let adapter = b"TTTTCCACGGGGATACTACTG";
        let mut r = ReadRec::from_seq(seq);
        let hit = trim_by_sequence(&mut r, &PreparedAdapter::new(adapter), 4, &mut AdapterScratch::default());
        assert!(hit.is_some());
        assert_eq!(r.bases(), b"TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAA");
    }

    /// fastp's own `AdapterTrimmer::test()`, multi-sequence case.
    #[test]
    fn fastp_trim_by_multi_sequences_case() {
        let seq = b"TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAATTTTCCCCGGGGAAATTTCCCGGGAAATTTCCCGGGATCGATCGATCGATCGAATTCC";
        let adapters: Vec<PreparedAdapter> = [
            b"GCTAGCTAGCTAGCTA".as_slice(),
            b"AAATTTCCCGGGAAATTTCCCGGG".as_slice(),
            b"ATCGATCGATCGATCG".as_slice(),
            b"AATTCCGGAATTCCGG".as_slice(),
        ]
        .iter()
        .map(|a| PreparedAdapter::new(a))
        .collect();
        let mut r = ReadRec::from_seq(seq);
        trim_by_multi_sequences(&mut r, &adapters, &mut AdapterScratch::default());
        assert_eq!(
            r.bases(),
            b"TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAATTTTCCCCGGGG"
        );
    }

    #[test]
    fn trims_a_real_truseq_tail() {
        // 40 bp of insert then the start of TruSeq R1
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
        let adapter = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
        let mut r = ReadRec::from_seq(seq);
        let hit = trim_by_sequence(&mut r, &PreparedAdapter::new(adapter), 4, &mut AdapterScratch::default()).expect("adapter should be found");
        assert_eq!(r.bases(), b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
        assert_eq!(hit.removed, 33);
    }

    #[test]
    fn leaves_a_clean_read_alone() {
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let adapter = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
        let mut r = ReadRec::from_seq(seq);
        assert!(trim_by_sequence(&mut r, &PreparedAdapter::new(adapter), 4, &mut AdapterScratch::default()).is_none());
        assert_eq!(r.bases(), seq);
    }

    #[test]
    fn an_adapter_only_read_is_emptied() {
        let adapter = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
        let mut r = ReadRec::from_seq(adapter);
        let hit = trim_by_sequence(&mut r, &PreparedAdapter::new(adapter), 4, &mut AdapterScratch::default()).expect("should match at zero");
        assert!(r.is_empty(), "left {:?}", r.bases());
        assert!(hit.removed > 0);
    }

    #[test]
    fn match_req_scales_with_the_list() {
        assert_eq!(match_req_for(1), 4);
        assert_eq!(match_req_for(17), 5);
        assert_eq!(match_req_for(257), 6);
    }
}
