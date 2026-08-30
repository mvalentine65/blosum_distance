//! Unit tests for [`translate`]. Compiled as a child module, so
//! private items stay reachable through `use super::*`.

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
