use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

// NCBI codon tables from https://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt
// The 64-char order is T/C/A/G for each codon position.
const TABLE_1: &[u8; 64] = b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
const TABLE_2: &[u8; 64] = b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG";
const TABLE_3: &[u8; 64] = b"FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
const TABLE_4: &[u8; 64] = b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
const TABLE_5: &[u8; 64] = b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG";
const TABLE_6: &[u8; 64] = b"FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
const TABLE_9: &[u8; 64] = b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG";
const TABLE_10: &[u8; 64] = b"FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
const TABLE_11: &[u8; 64] = b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
const TABLE_12: &[u8; 64] = b"FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
const TABLE_13: &[u8; 64] = b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG";
const TABLE_14: &[u8; 64] = b"FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG";
const TABLE_15: &[u8; 64] = b"FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
const TABLE_16: &[u8; 64] = b"FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
const TABLE_21: &[u8; 64] = b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG";
const TABLE_22: &[u8; 64] = b"FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
const TABLE_23: &[u8; 64] = b"FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
const TABLE_24: &[u8; 64] = b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG";
const TABLE_25: &[u8; 64] = b"FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
const TABLE_26: &[u8; 64] = b"FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
const TABLE_27: &[u8; 64] = b"FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
const TABLE_28: &[u8; 64] = b"FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
const TABLE_29: &[u8; 64] = b"FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
const TABLE_30: &[u8; 64] = b"FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
const TABLE_31: &[u8; 64] = b"FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
const TABLE_32: &[u8; 64] = b"FFLLSSSSYY*WCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
const TABLE_33: &[u8; 64] = b"FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG";

fn genetic_code_table(table: u8) -> Option<&'static [u8; 64]> {
    match table {
        1 => Some(TABLE_1),
        2 => Some(TABLE_2),
        3 => Some(TABLE_3),
        4 => Some(TABLE_4),
        5 => Some(TABLE_5),
        6 => Some(TABLE_6),
        9 => Some(TABLE_9),
        10 => Some(TABLE_10),
        11 => Some(TABLE_11),
        12 => Some(TABLE_12),
        13 => Some(TABLE_13),
        14 => Some(TABLE_14),
        15 => Some(TABLE_15),
        16 => Some(TABLE_16),
        21 => Some(TABLE_21),
        22 => Some(TABLE_22),
        23 => Some(TABLE_23),
        24 => Some(TABLE_24),
        25 => Some(TABLE_25),
        26 => Some(TABLE_26),
        27 => Some(TABLE_27),
        28 => Some(TABLE_28),
        29 => Some(TABLE_29),
        30 => Some(TABLE_30),
        31 => Some(TABLE_31),
        32 => Some(TABLE_32),
        33 => Some(TABLE_33),
        _ => None,
    }
}

#[inline]
fn nucleotide_index_tcag(base: u8) -> Option<usize> {
    match base {
        b'T' | b't' | b'U' | b'u' => Some(0),
        b'C' | b'c' => Some(1),
        b'A' | b'a' => Some(2),
        b'G' | b'g' => Some(3),
        _ => None,
    }
}

#[inline]
fn translate_codon(codon: &[u8], table: &[u8; 64]) -> char {
    if let (Some(base1), Some(base2), Some(base3)) = (
        nucleotide_index_tcag(codon[0]),
        nucleotide_index_tcag(codon[1]),
        nucleotide_index_tcag(codon[2]),
    ) {
        table[(base1 * 16) + (base2 * 4) + base3] as char
    } else {
        'X'
    }
}

/// translate(sequence: str, table: int = 1) -> str
///
/// Translate a nucleotide sequence with an NCBI codon table id, similar to
/// Bio.Seq.translate(table=<id>) default behavior.
///
/// Notes:
/// - Stops are returned as `*`.
/// - `U/u` is treated as `T/t`.
/// - Ambiguous/invalid codons return `X`.
/// - Trailing partial codons are ignored.
#[pyfunction]
#[pyo3(signature = (sequence, table = None), text_signature = "(sequence, table=1)")]
pub fn translate(sequence: &str, table: Option<u8>) -> PyResult<String> {
    let table_id = table.unwrap_or(1);
    let table = genetic_code_table(table_id).ok_or_else(|| {
        PyValueError::new_err(format!("Unknown NCBI codon table id: {table_id}"))
    })?;

    let mut output = String::with_capacity(sequence.len() / 3);
    for codon in sequence.as_bytes().chunks_exact(3) {
        output.push(translate_codon(codon, table));
    }
    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::translate;

    #[test]
    fn translates_standard_table() {
        let seq = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG";
        let translated = translate(seq, Some(1)).unwrap();
        assert_eq!(translated, "MAIVMGR*KGAR*");
    }

    #[test]
    fn table_argument_changes_translation() {
        let seq = "ATAAGAAGATGA";
        assert_eq!(translate(seq, Some(1)).unwrap(), "IRR*");
        assert_eq!(translate(seq, Some(2)).unwrap(), "M**W");
    }

    #[test]
    fn lower_case_and_u_are_supported() {
        assert_eq!(translate("augGcc", Some(1)).unwrap(), "MA");
    }

    #[test]
    fn ambiguous_codon_returns_x() {
        assert_eq!(translate("ATNCCC", Some(1)).unwrap(), "XP");
    }

    #[test]
    fn trailing_partial_codon_is_ignored() {
        assert_eq!(translate("ATGCC", Some(1)).unwrap(), "M");
    }

    #[test]
    fn default_table_works_when_not_provided() {
        assert_eq!(translate("ATG", None).unwrap(), "M");
    }

    #[test]
    fn unknown_table_returns_error() {
        assert!(translate("ATG", Some(8)).is_err());
    }
}
