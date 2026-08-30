//! Unit tests for [`crate`]. Compiled as a child module, so
//! private items stay reachable through `use super::*`.

use crate::blosum62_distance;

#[test]
fn test_blosum62_gap_penalty() {
    let result = blosum62_distance(String::from("A"), String::from("-"));
    assert_eq!(result, 2.0)
}
