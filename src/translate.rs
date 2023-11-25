/// Maps an ASCII character to array index
///
/// A = 65, a = 97  => 0
/// C = 67, c = 99  => 1
/// G = 71, g = 103 => 2
/// T = 84, t = 116 => 3
/// U = 85, u = 117 => 3
const ASCII_TO_INDEX: [usize;128] = [
        4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4, // 0-15
        4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4, // 16-31
        4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4, // 32-47
        4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4, // 48-63
        4,0,4,1, 4,4,4,2, 4,4,4,4, 4,4,4,4, // 64-79 (65 = A, 67 = C, 71 = G)
        4,4,4,4, 3,4,4,4, 4,4,4,4, 4,4,4,4, // 80-95 (84 = T)
        4,0,4,1, 4,4,4,2, 4,4,4,4, 4,4,4,4, // 96-111 (97 = a, 99 = c, 103 = g)
        4,4,4,4, 3,4,4,4, 4,4,4,4, 4,4,4,4, // 112-127 (116 = t)
                                    ];

 const AA_TABLE: [[[char; 4]; 4]; 4] = [
    [
        ['K', 'N', 'K', 'N'], // AAA, AAC, AAG, AAT
        ['T', 'T', 'T', 'T'], // ACA, ACC, ACG, ACT
        ['R', 'S', 'R', 'S'], // AGA, AGC, AGG, AGT
        ['I', 'I', 'M', 'I'], // ATA, ATC, ATG, ATT
    ],
    [
        ['Q', 'H', 'Q', 'H'], // CAA, CAC, CAG, CAT
        ['P', 'P', 'P', 'P'], // CCA, CCC, CCG, CCT
        ['R', 'R', 'R', 'R'], // CGA, CGC, CGG, CGT
        ['L', 'L', 'L', 'L'], // CTA, CTC, CTG, CTT
    ],
    [
        ['E', 'D', 'E', 'D'], // GAA, GAC, GAG, GAT
        ['A', 'A', 'A', 'A'], // GCA, GCC, GCG, GCT
        ['G', 'G', 'G', 'G'], // GGA, GGC, GGG, GGT
        ['V', 'V', 'V', 'V'], // GTA, GTC, GTG, GTT
    ],
    [
        ['*', 'Y', '*', 'Y'], // TAA, TAC, TAG, TAT
        ['S', 'S', 'S', 'S'], // TCA, TCC, TCG, TCT
        ['*', 'C', 'W', 'C'], // TGA, TGC, TGG, TGT
        ['L', 'F', 'L', 'F'], // TTA, TTC, TTG, TTT
    ],
];
use pyo3::prelude::*;
/// translate(sequence: str) -> str
/// Takes an nt string and converts it to an aa string.
/// Outputs "-" instead of "*".
/// The input must have a length that is divisible by 3
/// and must only contain 'A','T','C', and 'G' as ascii characters.
/// If the input length is not divisible by 3, returns an empty string.
/// Panics if an invalid character is found.
#[pyfunction]
pub fn translate(sequence: String) -> String {
    let sequence = sequence.as_bytes();
    if sequence.len() %3 != 0 || sequence.len() == 0 {
        return String::from("");
    }
    let mut output = Vec::with_capacity(sequence.len());
    for i in (0..sequence.len()).step_by(3) {
        output.push(AA_TABLE[ASCII_TO_INDEX[sequence[i] as usize]][ASCII_TO_INDEX[sequence[i+1] as usize]][ASCII_TO_INDEX[sequence[i+2] as usize]]);
    }
    output.iter().collect()
}