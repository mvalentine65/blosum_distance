//! Profile HMMs read from BATH/HMMER3 text files (bathbuild output).
//!
//! Probabilities are stored as the file gives them (`-ln p`, `*` for p = 0),
//! converted to f32 the way HMMER's reader does (expf). MAP gives each match node's
//! column in the alignment the model was built from.

use std::collections::HashMap;

pub const K: usize = 20; // residues, HMMER order ACDEFGHIKLMNPQRSTVWY

/// Transition indices, in file order.
pub const MM: usize = 0;
pub const MI: usize = 1;
pub const MD: usize = 2;
pub const IM: usize = 3;
pub const II: usize = 4;
pub const DM: usize = 5;
pub const DD: usize = 6;

pub struct Hmm {
    pub name: String,
    pub m: usize,
    /// Node 0..=m; node 0 holds the begin node's insert emissions and transitions.
    pub mat: Vec<[f32; K]>,
    /// Insert emissions: the distinct rows and each node's row (most nodes share one).
    ins_rows: Vec<[f32; K]>,
    ins_at: Vec<u32>,
    pub t: Vec<[f32; 7]>,
    /// 1-based alignment column of each match node (0 if the file has none).
    pub map: Vec<usize>,
}

impl Hmm {
    pub fn ins(&self, k: usize) -> &[f32; K] { &self.ins_rows[self.ins_at[k] as usize] }
}

fn prob(tok: &str) -> Result<f32, String> {
    if tok == "*" {
        return Ok(0.0);
    }
    let v: f64 = tok.parse().map_err(|_| format!("bad probability field '{tok}'"))?;
    // as HMMER reads it: expf(-1.0 * atof(tok)), single precision
    Ok(((-v) as f32).exp())
}

fn fields<const N: usize>(line: &str) -> Result<[f32; N], String> {
    let mut out = [0f32; N];
    let mut it = line.split_whitespace();
    for (i, o) in out.iter_mut().enumerate() {
        *o = prob(it.next().ok_or_else(|| format!("expected {N} fields, got {i}"))?)?;
    }
    Ok(out)
}

/// Every model in the file, in file order. Chunks of whole records parse on
/// `threads` threads.
pub fn read_hmms(path: &str, threads: usize) -> Result<Vec<Hmm>, String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("{path}: {e}"))?;
    let n = threads.max(1);
    let mut cuts = vec![0usize];
    for c in 1..n {
        let want = (text.len() * c / n).max(*cuts.last().unwrap());
        match text[want..].find("\n//") {
            Some(p) => {
                let at = text[want + p + 1..].find('\n').map_or(text.len(), |q| want + p + 1 + q + 1);
                if at > *cuts.last().unwrap() { cuts.push(at); }
            }
            None => break,
        }
    }
    cuts.push(text.len());
    cuts.dedup();
    let parts: Vec<Result<Vec<Hmm>, String>> = std::thread::scope(|sc| {
        let hs: Vec<_> = cuts.windows(2).map(|w| { let t = &text[w[0]..w[1]]; sc.spawn(move || parse_hmms(t)) }).collect();
        hs.into_iter().map(|h| h.join().expect("model parser panicked")).collect()
    });
    let mut out = Vec::new();
    for p in parts { out.extend(p?); }
    Ok(out)
}

fn parse_hmms(text: &str) -> Result<Vec<Hmm>, String> {
    let lines: Vec<&str> = text.lines().collect();
    let mut out = Vec::new();
    let mut i = 0;
    while i < lines.len() {
        // header
        let (mut name, mut m) = (String::new(), 0usize);
        while i < lines.len() && !lines[i].starts_with("HMM ") {
            let l = lines[i];
            if let Some(v) = l.strip_prefix("NAME") {
                name = v.trim().to_string();
            } else if let Some(v) = l.strip_prefix("LENG") {
                m = v.trim().parse().map_err(|_| format!("bad LENG line: {l}"))?;
            }
            i += 1;
        }
        if i >= lines.len() {
            break;
        }
        i += 2; // "HMM" line and the transition labels
        if lines.get(i).map_or(false, |l| l.trim_start().starts_with("COMPO")) {
            i += 1;
        }
        let mut h = Hmm {
            name,
            m,
            mat: vec![[0f32; K]; m + 1],
            ins_rows: Vec::new(),
            ins_at: vec![0; m + 1],
            t: vec![[0f32; 7]; m + 1],
            map: vec![0; m + 1],
        };
        let mut seen: HashMap<[u32; K], u32> = HashMap::new();
        let mut ins_row = |h: &mut Hmm, row: [f32; K]| -> u32 {
            *seen.entry(row.map(f32::to_bits)).or_insert_with(|| { h.ins_rows.push(row); h.ins_rows.len() as u32 - 1 })
        };
        h.ins_at[0] = ins_row(&mut h, fields::<K>(lines[i])?);
        h.t[0] = fields::<7>(lines[i + 1])?;
        i += 2;
        for k in 1..=m {
            // node number, K match emissions, then MAP and annotation fields
            let ml = lines[i].trim_start();
            let rest = ml.find(char::is_whitespace).map_or("", |p| &ml[p..]);
            h.mat[k] = fields::<K>(rest)?;
            h.map[k] = rest.split_whitespace().nth(K).and_then(|s| s.parse().ok()).unwrap_or(0);
            h.ins_at[k] = ins_row(&mut h, fields::<K>(lines[i + 1])?);
            h.t[k] = fields::<7>(lines[i + 2])?;
            i += 3;
        }
        while i < lines.len() && !lines[i].starts_with("//") {
            i += 1;
        }
        i += 1;
        out.push(h);
    }
    Ok(out)
}
