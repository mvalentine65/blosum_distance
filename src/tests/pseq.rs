//! Unit tests for [`pseq`]. Compiled as a child module, so
//! private items stay reachable through `use super::*`.

use super::*;

fn round_trip(chunk: &[u8]) -> Vec<u8> {
    let stored = pack(chunk).unwrap();
    let mut out = Vec::new();
    unpack_into(&stored, false, &mut out).unwrap();
    assert_eq!(out, chunk);
    let mut up = Vec::new();
    unpack_into(&stored, true, &mut up).unwrap();
    assert_eq!(up, chunk.to_ascii_uppercase());
    stored
}

fn random_chunk(len: usize, seed: &mut u64, alphabet: &[u8]) -> Vec<u8> {
    (0..len)
        .map(|_| {
            *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            alphabet[(*seed >> 33) as usize % alphabet.len()]
        })
        .collect()
}

#[test]
fn round_trips_at_every_length_edge() {
    let mut seed = 7;
    for len in [1, 2, 3, 4, 5, 63, 64, 65, MAX_CHUNK - 1, MAX_CHUNK] {
        round_trip(&random_chunk(len, &mut seed, b"ACGT"));
    }
}

#[test]
fn acgt_only_chunk_is_a_quarter_size() {
    let stored = round_trip(&vec![b'G'; MAX_CHUNK]);
    assert_eq!(stored.len(), HEADER + MAX_CHUNK / 4);
}

#[test]
fn soft_mask_and_ambiguity_codes_survive() {
    let mut seed = 11;
    // mixed case, N runs, every IUPAC code and U in both cases
    round_trip(&random_chunk(5000, &mut seed, b"ACGTacgtNnRYSWKMBDHVUryswkmbdhvu"));
    round_trip(&random_chunk(MAX_CHUNK, &mut seed, b"AAAACCCCGGGGTTTTaaaaccccggggttttN"));
}

#[test]
fn whole_chunk_runs_use_the_65536_length() {
    let lower = round_trip(&vec![b'a'; MAX_CHUNK]);
    assert_eq!(lower.len(), HEADER + MAX_CHUNK / 4 + CASE_RUN);
    let ns = round_trip(&vec![b'N'; MAX_CHUNK]);
    assert_eq!(ns.len(), HEADER + MAX_CHUNK / 4 + EX_RUN);
    round_trip(&vec![b'n'; MAX_CHUNK]);
}

#[test]
fn runs_at_both_ends() {
    round_trip(b"NNacgtACGTnn");
    round_trip(b"aCGTRRt");
}

#[test]
fn pack_refuses_what_it_cannot_round_trip() {
    assert!(pack(b"").is_err());
    assert!(pack(&vec![b'A'; MAX_CHUNK + 1]).is_err());
    assert!(pack(b"ACG-T").is_err());
    assert!(pack(b"ACG\nT").is_err());
}

#[test]
fn unpack_refuses_malformed_values() {
    let good = pack(b"ACGTNNacgt").unwrap();
    let mut out = Vec::new();
    assert!(unpack_into(&good[..HEADER - 1], false, &mut out).is_err());
    assert!(unpack_into(&good[..good.len() - 1], false, &mut out).is_err());
    let mut longer = good.clone();
    longer.push(0);
    assert!(unpack_into(&longer, false, &mut out).is_err());
    // a run reaching past the chunk's end
    let mut bad = good.clone();
    let ex = HEADER + 3;
    bad[ex] = 9;
    assert!(unpack_into(&bad, false, &mut out).is_err());
    // zero bases
    let mut zero = good.clone();
    zero[..4].copy_from_slice(&0u32.to_le_bytes());
    assert!(unpack_into(&zero, false, &mut out).is_err());
}

#[test]
fn unpack_appends_after_earlier_chunks() {
    let (a, b) = (pack(b"ACGTa").unwrap(), pack(b"nNTT").unwrap());
    let mut out = Vec::new();
    unpack_into(&a, false, &mut out).unwrap();
    unpack_into(&b, false, &mut out).unwrap();
    assert_eq!(out, b"ACGTanNTT");
}

#[test]
fn runs_at_every_offset_round_trip() {
    // random ACGT with lowercase and non-ACGT runs of random lengths starting anywhere
    let mut seed = 99u64;
    let mut rand = move |m: usize| {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (seed >> 33) as usize % m
    };
    for _ in 0..300 {
        let len = 1 + rand(3000);
        let mut chunk: Vec<u8> = (0..len).map(|_| b"ACGT"[rand(4)]).collect();
        for _ in 0..rand(12) {
            let s = rand(len);
            let e = (s + 1 + rand(40)).min(len);
            let fill = b"NNNRYU"[rand(6)];
            for c in &mut chunk[s..e] {
                *c = fill;
            }
        }
        for _ in 0..rand(20) {
            let s = rand(len);
            let e = (s + 1 + rand(60)).min(len);
            for c in &mut chunk[s..e] {
                *c |= 0x20;
            }
        }
        round_trip(&chunk);
    }
}

#[test]
fn word_path_packs_like_the_byte_path() {
    for x in [*b"ACGTACGT", *b"acgtACGT", *b"TTTTGGGG", *b"gcatgcat"] {
        let want: Vec<u8> = x.chunks(4).map(|f| f.iter().enumerate().fold(0u8, |a, (j, &b)| a | (CLASS[b as usize] & 3) << (2 * j))).collect();
        assert_eq!(pack8(u64::from_le_bytes(x)).unwrap().to_le_bytes().to_vec(), want);
    }
    assert!(pack8(u64::from_le_bytes(*b"ACGTNCGT")).is_none());
    assert!(pack8(u64::from_le_bytes(*b"ACGTACG\0")).is_none());
}
