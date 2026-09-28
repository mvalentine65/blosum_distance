//! GFF3 of the final gene copies: one gene feature per chain with exons, then
//! its exons in gene order, with refined splice sites where refine moved them.

use super::chain::{Chain, SRC_NAME};
use super::pseudo::Pseudo;
use super::refine::Refine;
use super::score::ExonScore;
use std::fmt::Write as _;

fn esc(s: &str) -> String {
    let mut o = String::with_capacity(s.len());
    for c in s.chars() {
        if matches!(c, ';' | '=' | '&' | ',' | '%' | '\t') { let _ = write!(o, "%{:02X}", c as u32); } else { o.push(c); }
    }
    o
}

pub fn write_gff(cs: &[Chain], rf: Option<&Refine>, ps: Option<&[Pseudo]>, es: Option<&[Vec<ExonScore>]>, min_strict: i64) -> String {
    let mut out = String::from("##gff-version 3\n");
    for (i, c) in cs.iter().enumerate() {
        if c.ex.is_empty() { continue; }
        let spans: Vec<(i64, i64)> = (0..c.ex.len())
            .map(|j| match rf { Some(r) => r.span(c, i, j), None => (c.ex[j].start, c.ex[j].end) })
            .collect();
        let lo = spans.iter().map(|x| x.0).min().unwrap();
        let hi = spans.iter().map(|x| x.1).max().unwrap();
        let _ = write!(out, "{}\texonfill\tgene\t{}\t{}\t.\t{}\t.\tID={};Name={};nodes={}-{}", c.scaffold, lo, hi,
                       c.strand as char, esc(&c.gene), esc(&c.model), c.klo, c.khi);
        if let Some(ps) = ps {
            let p = &ps[i];
            let _ = write!(out, ";pseudogene={}", if c.ex.len() > 1 && p.called(min_strict) { "yes" } else { "no" });
            if p.fs_tight != 0 || p.stop_tight != 0 || p.nnonc != 0 {
                let _ = write!(out, ";disablements=fs:{},stop:{},noncanonical:{}", p.fs_tight, p.stop_tight, p.nnonc);
            }
        }
        out.push('\n');
        for (j, e) in c.ex.iter().enumerate() {
            let _ = write!(out, "{}\texonfill\texon\t{}\t{}\t.\t{}\t.\tID={}.e{};Parent={};source={};nodes={}-{}", c.scaffold,
                           spans[j].0, spans[j].1, c.strand as char, esc(&c.gene), j + 1, esc(&c.gene), SRC_NAME[e.src as usize], e.k1, e.k2);
            if let Some(r) = rf {
                let x = &r.ex[i][j];
                if x.acc_g != 0 { let _ = write!(out, ";acceptor={}", String::from_utf8_lossy(&x.acc_site)); }
                if x.don_g != 0 { let _ = write!(out, ";donor={}", String::from_utf8_lossy(&x.don_site)); }
            }
            if let Some(x) = es.map(|v| &v[i][j]) {
                if x.bits.is_finite() { let _ = write!(out, ";bits={:.1};rev={:.1}", x.bits, x.rev); }
            }
            out.push('\n');
        }
    }
    out
}
