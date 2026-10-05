//! Codon alignment of nucleotide sequences threaded onto their aligned peptides.
//! Moved from the pr2codon package (mvalentine65/pn2codon af0f3b4, 1.1.22).
//!
//! pn2codon began as a port of the pn2codon routine of PAL2NAL (pal2nal.pl,
//! Mikita Suyama) and has since been rewritten: it checks consecutive codons
//! against the NCBI tables instead of searching with PAL2NAL's
//! reverse-translation patterns. Suyama, Torrents & Bork (2006), Nucleic Acids
//! Research 34:W609-W612.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::cell::{Cell, RefCell};
use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::sync::Arc;

type GeneTable = HashMap<char, Vec<String>>;

const VALID_PEPS: &[char] = &[
    'A', 'L', 'W', 'Q', 'Y', 'E', 'C', 'D', 'F', 'G', 'H', 'I', 'M', 'K', 'P', 'R', 'S', 'V',
    'N', 'T', '*', '-', 'B', 'J', 'Z', 'X',
];

fn table_signature(table_num: i32) -> Option<&'static str> {
    let table = crate::translate::genetic_code_table(u8::try_from(table_num).ok()?)?;
    std::str::from_utf8(table).ok()
}

fn normalize_base(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'U' => b'T',
        normalized => normalized,
    }
}

fn iupac_expansions(base: u8) -> &'static [u8] {
    match normalize_base(base) {
        b'R' => b"AG",
        b'Y' => b"CT",
        b'S' => b"GC",
        b'W' => b"AT",
        b'K' => b"GT",
        b'M' => b"AC",
        b'B' => b"CGT",
        b'D' => b"AGT",
        b'H' => b"ACT",
        b'V' => b"ACG",
        b'N' => b"ACGT",
        _ => b"",
    }
}

fn base_rank(base: u8) -> Option<usize> {
    match normalize_base(base) {
        b'T' => Some(0),
        b'C' => Some(1),
        b'A' => Some(2),
        b'G' => Some(3),
        _ => None,
    }
}

fn codon_index(codon: &[u8]) -> Option<usize> {
    if codon.len() != 3 {
        return None;
    }
    Some(base_rank(codon[0])? * 16 + base_rank(codon[1])? * 4 + base_rank(codon[2])?)
}

fn aa_matches_table_symbol(table_aa: u8, aa: char) -> bool {
    match aa {
        'B' => table_aa == b'N' || table_aa == b'D',
        'J' => table_aa == b'L' || table_aa == b'I',
        'Z' => table_aa == b'Q' || table_aa == b'E',
        _ => table_aa == aa as u8,
    }
}

fn codon_matches_signature(signature: &str, aa: char, codon: &str) -> bool {
    match codon_index(codon.as_bytes()) {
        Some(idx) => match signature.as_bytes().get(idx) {
            Some(table_aa) => aa_matches_table_symbol(*table_aa, aa),
            None => false,
        },
        None => false,
    }
}

fn fill_possible_bases(base: u8, out: &mut [u8; 4]) -> usize {
    match normalize_base(base) {
        b'A' => {
            out[0] = b'A';
            1
        }
        b'C' => {
            out[0] = b'C';
            1
        }
        b'G' => {
            out[0] = b'G';
            1
        }
        b'T' => {
            out[0] = b'T';
            1
        }
        b'R' => {
            out[0] = b'A';
            out[1] = b'G';
            2
        }
        b'Y' => {
            out[0] = b'C';
            out[1] = b'T';
            2
        }
        b'S' => {
            out[0] = b'G';
            out[1] = b'C';
            2
        }
        b'W' => {
            out[0] = b'A';
            out[1] = b'T';
            2
        }
        b'K' => {
            out[0] = b'G';
            out[1] = b'T';
            2
        }
        b'M' => {
            out[0] = b'A';
            out[1] = b'C';
            2
        }
        b'B' => {
            out[0] = b'C';
            out[1] = b'G';
            out[2] = b'T';
            3
        }
        b'D' => {
            out[0] = b'A';
            out[1] = b'G';
            out[2] = b'T';
            3
        }
        b'H' => {
            out[0] = b'A';
            out[1] = b'C';
            out[2] = b'T';
            3
        }
        b'V' => {
            out[0] = b'A';
            out[1] = b'C';
            out[2] = b'G';
            3
        }
        b'N' => {
            out[0] = b'A';
            out[1] = b'C';
            out[2] = b'G';
            out[3] = b'T';
            4
        }
        _ => 0,
    }
}

fn ambiguous_triplet_matches_signature(signature: &str, aa: char, triplet: &str) -> bool {
    let bytes = triplet.as_bytes();
    if bytes.len() != 3 {
        return false;
    }

    let signature_bytes = signature.as_bytes();
    let mut p0 = [0_u8; 4];
    let mut p1 = [0_u8; 4];
    let mut p2 = [0_u8; 4];
    let l0 = fill_possible_bases(bytes[0], &mut p0);
    let l1 = fill_possible_bases(bytes[1], &mut p1);
    let l2 = fill_possible_bases(bytes[2], &mut p2);

    if l0 == 0 || l1 == 0 || l2 == 0 {
        return false;
    }

    for b0 in p0.iter().take(l0) {
        for b1 in p1.iter().take(l1) {
            for b2 in p2.iter().take(l2) {
                let idx = match codon_index(&[*b0, *b1, *b2]) {
                    Some(i) => i,
                    None => continue,
                };
                if let Some(table_aa) = signature_bytes.get(idx) {
                    if aa_matches_table_symbol(*table_aa, aa) {
                        return true;
                    }
                }
            }
        }
    }

    false
}

fn recurse(triplet: &[u8], working: &mut [u8], index: usize, output: &mut HashSet<String>) {
    if index == triplet.len() {
        output.insert(String::from_utf8(working.to_vec()).unwrap_or_default());
        return;
    }

    let this_char = normalize_base(triplet[index]);
    working[index] = this_char;
    recurse(triplet, working, index + 1, output);

    for replacement in iupac_expansions(this_char) {
        working[index] = *replacement;
        recurse(triplet, working, index + 1, output);
    }
}

pub fn make_iupac_set(triplet: &[u8]) -> HashSet<String> {
    let mut output = HashSet::new();
    let mut working = vec![0_u8; triplet.len()];
    recurse(triplet, &mut working, 0, &mut output);
    output
}

fn iupac_matches(pattern_base: u8, concrete_base: u8) -> bool {
    let pattern = normalize_base(pattern_base);
    let concrete = normalize_base(concrete_base);
    match pattern {
        b'A' | b'C' | b'G' | b'T' => concrete == pattern,
        b'R' => concrete == b'A' || concrete == b'G',
        b'Y' => concrete == b'C' || concrete == b'T',
        b'S' => concrete == b'G' || concrete == b'C',
        b'W' => concrete == b'A' || concrete == b'T',
        b'K' => concrete == b'G' || concrete == b'T',
        b'M' => concrete == b'A' || concrete == b'C',
        b'B' => concrete == b'C' || concrete == b'G' || concrete == b'T',
        b'D' => concrete == b'A' || concrete == b'G' || concrete == b'T',
        b'H' => concrete == b'A' || concrete == b'C' || concrete == b'T',
        b'V' => concrete == b'A' || concrete == b'C' || concrete == b'G',
        b'N' => matches!(concrete, b'A' | b'C' | b'G' | b'T'),
        _ => false,
    }
}

fn has_iupac_match<'a, I>(original_triplet: &str, taxa: I) -> bool
where
    I: IntoIterator<Item = &'a str>,
{
    let pattern = original_triplet.as_bytes();
    if pattern.len() != 3 {
        return false;
    }

    for triplet in taxa {
        let candidate = triplet.as_bytes();
        if candidate.len() != 3 {
            continue;
        }

        if pattern
            .iter()
            .zip(candidate.iter())
            .all(|(pattern_base, concrete_base)| iupac_matches(*pattern_base, *concrete_base))
        {
            return true;
        }
    }
    false
}

#[pyfunction]
pub fn attempt_iupac_substitution(original_triplet: &str, taxa: Vec<String>) -> Option<String> {
    let possible_subs = make_iupac_set(original_triplet.as_bytes());
    for triplet in taxa {
        if possible_subs.contains(triplet.as_str()) {
            return Some(original_triplet.to_string());
        }
    }
    None
}

#[derive(Clone)]
struct AminoAcidTranslator {
    sequence_index: usize,
    aa_source_label: Arc<str>,
    aa_header: String,
    amino_acid: String,
    nt_source_label: Arc<str>,
    nt_header: String,
    nucleotide: String,
    has_reported_error: Cell<bool>,
    error_message: RefCell<Option<String>>,
}

fn truncate_header(header: &str) -> String {
    if header.len() > 150 {
        format!("{}...", &header[0..150])
    } else {
        header.to_string()
    }
}

fn source_label_from_path(path: &str, fallback: &str) -> String {
    let trimmed = path.trim();
    if trimmed.is_empty() {
        return fallback.to_string();
    }

    let file_name = Path::new(trimmed)
        .file_name()
        .and_then(|value| value.to_str())
        .filter(|value| !value.is_empty())
        .unwrap_or(trimmed);

    if file_name.to_ascii_lowercase().ends_with(".fa") {
        file_name.to_string()
    } else {
        format!("{}.fa", file_name)
    }
}

fn format_error_block(title: &str, details: &str) -> String {
    let mut out = String::new();
    out.push_str("========================================\n");
    out.push_str(&format!("pn2codon ERROR: {}\n", title));
    out.push_str("========================================\n");
    if !details.is_empty() {
        out.push_str(details);
        if !details.ends_with('\n') {
            out.push('\n');
        }
    }
    out.push_str("========================================\n");
    out
}

fn clamp_window(len: usize, center: usize, radius: usize) -> (usize, usize) {
    if len == 0 {
        return (0, 0);
    }
    let clamped_center = center.min(len - 1);
    let start = clamped_center.saturating_sub(radius);
    let end = (clamped_center + radius + 1).min(len);
    (start, end)
}

#[derive(Clone, Copy)]
enum NtTrackMode<'a> {
    None,
    Signature(&'a str),
    Table(&'a GeneTable),
}

fn spaced_aa_track(aas: &str) -> String {
    if aas.is_empty() {
        return String::new();
    }

    let mut out = String::new();
    out.push(' ');
    for (idx, aa) in aas.chars().enumerate() {
        if idx > 0 {
            out.push_str("  ");
        }
        out.push(aa);
    }
    out
}

fn translate_nt_window_with_signature(signature: &str, nt_window: &str) -> String {
    let mut out = String::new();
    for chunk in nt_window.as_bytes().chunks(3) {
        if chunk.len() != 3 {
            break;
        }

        let aa = codon_index(chunk)
            .and_then(|idx| signature.as_bytes().get(idx).copied())
            .map(char::from)
            .unwrap_or('X');
        out.push(aa);
    }
    out
}

fn translate_nt_window_with_table(gene_table: &GeneTable, nt_window: &str) -> String {
    let mut aa_keys: Vec<char> = gene_table.keys().copied().collect();
    aa_keys.sort_unstable();

    let mut out = String::new();
    for chunk in nt_window.as_bytes().chunks(3) {
        if chunk.len() != 3 {
            break;
        }

        let codon = match std::str::from_utf8(chunk) {
            Ok(value) => value,
            Err(_) => {
                out.push('X');
                continue;
            }
        };

        let mut mapped = 'X';
        for aa in aa_keys.iter() {
            if let Some(codons) = gene_table.get(aa) {
                if codons.iter().any(|triplet| triplet == codon)
                    || has_iupac_match(codon, codons.iter().map(String::as_str))
                {
                    mapped = *aa;
                    break;
                }
            }
        }
        out.push(mapped);
    }
    out
}

fn codon_matches_expected_in_mode(nt_track_mode: NtTrackMode<'_>, expected_aa: char, codon: &str) -> bool {
    match nt_track_mode {
        NtTrackMode::None => true,
        NtTrackMode::Signature(signature) => {
            codon_matches_signature(signature, expected_aa, codon)
                || ambiguous_triplet_matches_signature(signature, expected_aa, codon)
        }
        NtTrackMode::Table(gene_table) => match gene_table.get(&expected_aa) {
            Some(taxa) => {
                taxa.iter().any(|triplet| triplet == codon)
                    || has_iupac_match(codon, taxa.iter().map(String::as_str))
            }
            None => false,
        },
    }
}

fn residue_index_for_alignment_position(aligned_aa: &str, aa_alignment_index: usize) -> usize {
    let residue_count = aligned_aa
        .chars()
        .take(aa_alignment_index + 1)
        .filter(|c| *c != '-' && !c.is_ascii_digit())
        .count();
    residue_count.saturating_sub(1)
}

fn format_seq_inconsistency_details(
    aa_source_label: &str,
    aa_id: &str,
    aa_seq: &str,
    nt_source_label: &str,
    nt_id: &str,
    nt_seq: &str,
    aa_center_index: usize,
    nt_track_mode: NtTrackMode<'_>,
    mismatch_nt_base_index: Option<usize>,
) -> String {
    let aa_center = if aa_seq.is_empty() {
        0
    } else {
        aa_center_index.min(aa_seq.len().saturating_sub(1))
    };
    let (aa_start, aa_end) = clamp_window(aa_seq.len(), aa_center, 10);
    let mut nt_start = aa_start.saturating_mul(3);
    let nt_end = (aa_end.saturating_mul(3)).min(nt_seq.len());
    if nt_start > nt_end {
        nt_start = nt_end;
    }
    let aa_window = if aa_start < aa_end {
        &aa_seq[aa_start..aa_end]
    } else {
        ""
    };
    let nt_window = if nt_start < nt_end {
        &nt_seq[nt_start..nt_end]
    } else {
        ""
    };
    let aa_prefix = if aa_start > 0 { "." } else { "" };
    let aa_suffix = if aa_end < aa_seq.len() { "." } else { "" };
    let nt_prefix = if nt_start > 0 { "..." } else { "" };
    let nt_suffix = if nt_end < nt_seq.len() { "..." } else { "" };

    let mut out = String::new();
    out.push_str(&format!("Peptide header ({}) : {}\n", aa_source_label, aa_id));
    out.push_str(&format!(
        "Nucleotide header ({}) : {}\n\n",
        nt_source_label, nt_id
    ));
    out.push_str(&format!(
        ">{} {}-{}\n",
        aa_id,
        aa_start.saturating_add(1),
        aa_end
    ));
    if aa_window.is_empty() {
        out.push_str("(empty)\n");
    } else {
        out.push_str(aa_prefix);
        out.push_str(aa_window);
        out.push_str(aa_suffix);
        out.push('\n');
        let aa_track = spaced_aa_track(aa_window);
        if !aa_track.is_empty() {
            out.push_str(&aa_track);
            out.push('\n');
        }
    }
    out.push_str(&format!(
        ">{} {}-{}\n",
        nt_id,
        nt_start.saturating_add(1),
        nt_end
    ));
    if nt_window.is_empty() {
        out.push_str("(empty)\n");
    } else {
        out.push_str(nt_prefix);
        out.push_str(nt_window);
        out.push_str(nt_suffix);
        out.push('\n');

        let translated_nt_aas = match nt_track_mode {
            NtTrackMode::None => String::new(),
            NtTrackMode::Signature(signature) => translate_nt_window_with_signature(signature, nt_window),
            NtTrackMode::Table(gene_table) => translate_nt_window_with_table(gene_table, nt_window),
        };

        let nt_track = spaced_aa_track(&translated_nt_aas);
        if !nt_track.is_empty() {
            out.push_str(&" ".repeat(nt_prefix.len()));
            out.push_str(&nt_track);
            out.push('\n');
        }

        if mismatch_nt_base_index.is_some() && !matches!(nt_track_mode, NtTrackMode::None) {
            let expected_aas: Vec<char> = aa_window.chars().collect();
            let codon_count = expected_aas.len().min(nt_window.len() / 3);
            let mut markers = vec![' '; nt_window.len()];

            for codon_idx in 0..codon_count {
                let codon_start = codon_idx * 3;
                let codon = &nt_window[codon_start..codon_start + 3];
                let expected_aa = expected_aas[codon_idx];
                if !codon_matches_expected_in_mode(nt_track_mode, expected_aa, codon) {
                    markers[codon_start] = '_';
                    markers[codon_start + 1] = '_';
                    markers[codon_start + 2] = '_';
                }
            }

            if markers.iter().any(|c| *c == '_') {
                let marker_line: String = markers.into_iter().collect();
                out.push_str(&" ".repeat(nt_prefix.len()));
                out.push_str(&marker_line);
                out.push('\n');
            }
        }
    }
    out
}

fn supported_table_numbers() -> String {
    (1..=u8::MAX)
        .filter(|id| crate::translate::genetic_code_table(*id).is_some())
        .map(|id| id.to_string())
        .collect::<Vec<String>>()
        .join(", ")
}

impl AminoAcidTranslator {
    fn new(
        sequence_index: usize,
        aa_source_label: Arc<str>,
        aa_header: String,
        amino_acid: String,
        nt_source_label: Arc<str>,
        nt_header: String,
        nucleotide: String,
    ) -> Self {
        Self {
            sequence_index,
            aa_source_label,
            aa_header,
            amino_acid,
            nt_source_label,
            nt_header,
            nucleotide,
            has_reported_error: Cell::new(false),
            error_message: RefCell::new(None),
        }
    }

    fn do_checks(&self) {
        if self.aa_header != self.nt_header {
            let details = format!(
                "Sequence index: {}\nExpected header ({}): \"{}\"\nheader found ({}): \"{}\"",
                self.sequence_index,
                self.aa_source_label.as_ref(),
                truncate_header(&self.aa_header),
                self.nt_source_label.as_ref(),
                truncate_header(&self.nt_header)
            );
            self.report_error(format_error_block(
                "Header mismatch between peptide and nucleotide records.",
                &details,
            ));
        }

        let aa_triplet_len = self.amino_acid.chars().filter(|c| *c != '-').count() * 3;
        let nt_len = self.nucleotide.len();

        if nt_len != aa_triplet_len {
            let aa_seq: String = self
                .amino_acid
                .chars()
                .filter(|c| *c != '-' && !c.is_ascii_digit())
                .collect();
            let nt_seq = self.nucleotide.clone();
            let aa_center = if aa_seq.is_empty() {
                0
            } else {
                (nt_len.saturating_sub(1) / 3).min(aa_seq.len().saturating_sub(1))
            };
            let details = format_seq_inconsistency_details(
                self.aa_source_label.as_ref(),
                &truncate_header(&self.aa_header),
                &aa_seq,
                self.nt_source_label.as_ref(),
                &truncate_header(&self.nt_header),
                &nt_seq,
                aa_center,
                NtTrackMode::None,
                None,
            );
            self.report_error(format_error_block(
                "Peptide and nucleotide lengths are inconsistent.",
                &details,
            ));
        }
    }

    fn streamline(&mut self) {
        self.amino_acid = self
            .amino_acid
            .trim()
            .chars()
            .map(|c| {
                // '.' is an alignment gap (HMMER/needletail emit '.' for
                // insert-state gaps, '-' for match-state gaps). Normalise it to
                // '-' so the length check in `do_checks` — which counts
                // `filter(|c| *c != '-')` — excludes it, exactly as the
                // nucleotide track strips both '-' and '.'. Without this, '.'
                // falls through to the `else` branch and becomes 'X', which is
                // then miscounted as a residue and yields a spurious
                // "Peptide and nucleotide lengths are inconsistent." error.
                if c == '.' {
                    return '-';
                }
                let upper = c.to_ascii_uppercase();
                if VALID_PEPS.contains(&upper) {
                    upper
                } else {
                    'X'
                }
            })
            .collect();

        self.nucleotide = self
            .nucleotide
            .trim()
            .to_ascii_uppercase()
            .replace('-', "")
            .replace('.', "");
    }

    fn report_error(&self, message: String) {
        if self.has_reported_error.replace(true) {
            return;
        }
        self.error_message.replace(Some(message));
    }

    fn get_error_message(&self) -> Option<String> {
        self.error_message.borrow().clone()
    }

    fn error_out_mismatch(&self, aa_index: usize, nt_base_index: usize, nt_track_mode: NtTrackMode<'_>) {
        let aa_seq: String = self
            .amino_acid
            .chars()
            .filter(|c| *c != '-' && !c.is_ascii_digit())
            .collect();
        let nt_seq = self.nucleotide.clone();
        let aa_center = residue_index_for_alignment_position(&self.amino_acid, aa_index);
        let details = format_seq_inconsistency_details(
            self.aa_source_label.as_ref(),
            &truncate_header(&self.aa_header),
            &aa_seq,
            self.nt_source_label.as_ref(),
            &truncate_header(&self.nt_header),
            &nt_seq,
            aa_center,
            nt_track_mode,
            Some(nt_base_index),
        );
        self.report_error(format_error_block(
            "Peptide and nucleotide sequences are inconsistent.",
            &details,
        ));
    }

    fn reverse_translate_and_compare_with_table(&self, gene_table: &GeneTable) -> String {
        let mut compare_triplets = self.nucleotide.as_bytes().chunks(3);
        let mut nt_triplet_index = 0_usize;
        let mut output = String::with_capacity(self.nucleotide.len());

        for (aa_index, aa) in self.amino_acid.chars().enumerate() {
            if aa == '-' {
                output.push_str("---");
                continue;
            }

            if aa.is_ascii_digit() {
                if let Some(digit) = aa.to_digit(10) {
                    output.push_str(&".".repeat(digit as usize));
                }
                continue;
            }

            let taxa = match gene_table.get(&aa) {
                Some(codons) => codons,
                None => {
                    let mut supported: Vec<char> = gene_table.keys().copied().collect();
                    supported.sort_unstable();
                    let supported_list = supported
                        .iter()
                        .map(char::to_string)
                        .collect::<Vec<String>>()
                        .join(", ");
                    let details = format!(
                        "Amino acid         : '{}'\nAlignment position : {}\nValid symbols      : {}",
                        aa,
                        aa_index + 1,
                        supported_list
                    );
                    self.report_error(format_error_block(
                        "Amino acid is missing from custom codon table.",
                        &details,
                    ));
                    return String::new();
                }
            };

            let nt_base_index = nt_triplet_index * 3;
            let original_triplet = match compare_triplets.next() {
                Some(chunk) if chunk.len() == 3 => match std::str::from_utf8(chunk) {
                    Ok(triplet) => {
                        nt_triplet_index += 1;
                        triplet
                    }
                    Err(_) => {
                        self.error_out_mismatch(
                            aa_index,
                            nt_base_index,
                            NtTrackMode::Table(gene_table),
                        );
                        return String::new();
                    }
                },
                _ => {
                    self.error_out_mismatch(aa_index, nt_base_index, NtTrackMode::Table(gene_table));
                    return String::new();
                }
            };

            if original_triplet.contains('N') || aa == 'X' {
                output.push_str(original_triplet);
                continue;
            }

            if taxa.iter().any(|triplet| triplet == original_triplet)
                || has_iupac_match(original_triplet, taxa.iter().map(String::as_str))
            {
                output.push_str(original_triplet);
            } else {
                self.error_out_mismatch(aa_index, nt_base_index, NtTrackMode::Table(gene_table));
                return String::new();
            }
        }

        output
    }

    fn reverse_translate_and_compare_with_signature(&self, signature: &str) -> String {
        let mut compare_triplets = self.nucleotide.as_bytes().chunks(3);
        let mut nt_triplet_index = 0_usize;
        let mut output = String::with_capacity(self.nucleotide.len());

        for (aa_index, aa) in self.amino_acid.chars().enumerate() {
            if aa == '-' {
                output.push_str("---");
                continue;
            }

            if aa.is_ascii_digit() {
                if let Some(digit) = aa.to_digit(10) {
                    output.push_str(&".".repeat(digit as usize));
                }
                continue;
            }

            let nt_base_index = nt_triplet_index * 3;
            let original_triplet = match compare_triplets.next() {
                Some(chunk) if chunk.len() == 3 => match std::str::from_utf8(chunk) {
                    Ok(triplet) => {
                        nt_triplet_index += 1;
                        triplet
                    }
                    Err(_) => {
                        self.error_out_mismatch(
                            aa_index,
                            nt_base_index,
                            NtTrackMode::Signature(signature),
                        );
                        return String::new();
                    }
                },
                _ => {
                    self.error_out_mismatch(
                        aa_index,
                        nt_base_index,
                        NtTrackMode::Signature(signature),
                    );
                    return String::new();
                }
            };

            if original_triplet.contains('N') || aa == 'X' {
                output.push_str(original_triplet);
                continue;
            }

            if codon_matches_signature(signature, aa, original_triplet)
                || ambiguous_triplet_matches_signature(signature, aa, original_triplet)
            {
                output.push_str(original_triplet);
            } else {
                self.error_out_mismatch(
                    aa_index,
                    nt_base_index,
                    NtTrackMode::Signature(signature),
                );
                return String::new();
            }
        }

        output
    }
}

fn translate_record_with_table(
    gene_table: &GeneTable,
    sequence_index: usize,
    aa_source_label: &Arc<str>,
    nt_source_label: &Arc<str>,
    aa_header: String,
    aa: String,
    nt_header: String,
    nt: String,
) -> Result<String, String> {
    let mut translator = AminoAcidTranslator::new(
        sequence_index,
        Arc::clone(aa_source_label),
        aa_header,
        aa,
        Arc::clone(nt_source_label),
        nt_header,
        nt,
    );
    translator.streamline();
    translator.do_checks();
    if translator.has_reported_error.get() {
        return Err(translator.get_error_message().unwrap_or_else(|| {
            format_error_block("Peptide and nucleotide sequences are inconsistent.", "")
        }));
    }
    let codon = translator.reverse_translate_and_compare_with_table(gene_table);
    if translator.has_reported_error.get() {
        return Err(translator.get_error_message().unwrap_or_else(|| {
            format_error_block("Peptide and nucleotide sequences are inconsistent.", "")
        }));
    }
    Ok(codon)
}

fn translate_record_with_signature(
    signature: &str,
    sequence_index: usize,
    aa_source_label: &Arc<str>,
    nt_source_label: &Arc<str>,
    aa_header: String,
    aa: String,
    nt_header: String,
    nt: String,
) -> Result<String, String> {
    let mut translator = AminoAcidTranslator::new(
        sequence_index,
        Arc::clone(aa_source_label),
        aa_header,
        aa,
        Arc::clone(nt_source_label),
        nt_header,
        nt,
    );
    translator.streamline();
    translator.do_checks();
    if translator.has_reported_error.get() {
        return Err(translator.get_error_message().unwrap_or_else(|| {
            format_error_block("Peptide and nucleotide sequences are inconsistent.", "")
        }));
    }
    let codon = translator.reverse_translate_and_compare_with_signature(signature);
    if translator.has_reported_error.get() {
        return Err(translator.get_error_message().unwrap_or_else(|| {
            format_error_block("Peptide and nucleotide sequences are inconsistent.", "")
        }));
    }
    Ok(codon)
}

#[pyfunction]
pub fn pn2codon(
    _file_steem: String,
    aa_path: String,
    nt_path: String,
    table_num: i32,
    seqs: HashMap<String, ((String, String), (i32, String, String))>,
) -> PyResult<String> {
    let aa_source_label: Arc<str> = source_label_from_path(&aa_path, "aa.fa").into();
    let nt_source_label: Arc<str> = source_label_from_path(&nt_path, "nt.fa").into();

    let signature = match table_signature(table_num) {
        Some(sig) => sig,
        None => {
            let details = format!(
                "Requested table : {}\nSupported NCBI tables: {}",
                table_num,
                supported_table_numbers()
            );
            return Err(PyValueError::new_err(format_error_block(
                "Invalid codon table number.",
                &details,
            )));
        }
    };

    // In input order: a hash map's order, and with it the record an error names,
    // differs from call to call.
    let mut records: Vec<_> = seqs.into_iter().collect();
    records.sort_by(|a, b| a.1.1.0.cmp(&b.1.1.0).then_with(|| a.0.cmp(&b.0)));

    let mut file = String::new();
    for (index, (header, ((aa_header, aa), (_, nt_header, nt)))) in records.into_iter().enumerate() {
        let codon = translate_record_with_signature(
            signature,
            index + 1,
            &aa_source_label,
            &nt_source_label,
            aa_header,
            aa,
            nt_header,
            nt,
        )
        .map_err(PyValueError::new_err)?;
        file.push_str(&header);
        file.push('\n');
        file.push_str(&codon);
        file.push('\n');
    }
    Ok(file)
}

#[pyfunction]
pub fn pn2codon_original_args(
    _file_steem: String,
    aa_path: String,
    nt_path: String,
    gene_table: HashMap<char, Vec<String>>,
    seqs: HashMap<String, ((String, String), (String, String))>,
) -> PyResult<String> {
    let aa_source_label: Arc<str> = source_label_from_path(&aa_path, "aa.fa").into();
    let nt_source_label: Arc<str> = source_label_from_path(&nt_path, "nt.fa").into();

    // By header: these records carry no input position.
    let mut records: Vec<_> = seqs.into_iter().collect();
    records.sort_by(|a, b| a.0.cmp(&b.0));

    let mut file = String::new();
    for (index, (header, ((aa_header, aa), (nt_header, nt)))) in records.into_iter().enumerate() {
        let codon = translate_record_with_table(
            &gene_table,
            index + 1,
            &aa_source_label,
            &nt_source_label,
            aa_header,
            aa,
            nt_header,
            nt,
        )
        .map_err(PyValueError::new_err)?;
        file.push_str(&header);
        file.push('\n');
        file.push_str(&codon);
        file.push('\n');
    }
    Ok(file)
}
