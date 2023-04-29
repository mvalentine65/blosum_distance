use pyo3::prelude::*;
use serde::Deserialize;
use serde::Serialize;

#[pyclass]
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct ReferenceHit {
    target: String,
    sstart: i32,
    send: i32,
}

#[pymethods]
impl ReferenceHit {
    #[new]
    fn new(target: String, sstart: i32, send: i32) -> Self {
        Self {
            target,
            sstart,
            send,
        }
    }

    fn to_json(&self) -> String {
        serde_json::to_string(&self).unwrap()
    }
}

#[pyclass]
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct Hit {
    #[pyo3(get, set)]
    header: String,
    #[pyo3(get, set)]
    target: String,
    #[pyo3(get, set)]
    gene: Option<String>,
    #[pyo3(get, set)]
    frame: i8,
    #[pyo3(get, set)]
    evalue: f64,
    #[pyo3(get, set)]
    score: f32,
    #[pyo3(get, set)]
    qstart: i16,
    #[pyo3(get, set)]
    qend: i16,
    #[pyo3(get, set)]
    sstart: i32,
    #[pyo3(get, set)]
    send: i32,
    #[pyo3(get, set)]
    pident: f32,
    #[pyo3(get, set)]
    reftaxon: Option<String>,
    #[pyo3(get, set)]
    kick: bool,
    #[pyo3(get, set)]
    seq: Option<String>,
    #[pyo3(get, set)]
    length: i32,
    #[pyo3(get, set)]
    full_header: String,
    #[pyo3(get, set)]
    reference_hits: Vec<ReferenceHit>,
}

#[pymethods]
impl Hit {
    #[new]
    fn new(
        header: String,
        ref_header: String,
        frame: i8,
        evaule: f64,
        score: f32,
        qstart: i16,
        qend: i16,
        sstart: i32,
        send: i32,
        pident: f32,
    ) -> Self {
        let _qstart;
        let _qend;
        let _full_header;
        let frame = frame;
        let sstart = sstart;
        let send = send;
        if frame >= 0 {
            // normal values
            _qstart = qstart;
            _qend = qend;
            _full_header = format!("{}|[translate({})]", header, frame);
        } else {
            // frame < 0, so reverse and mark as revcomp
            _qstart = qend;
            _qend = qstart;
            _full_header = format!("{}|[revcomp]:[translate({})]", header, frame.abs());
        };
        Self {
            header: header,
            target: ref_header.clone(),
            gene: None,
            reftaxon: None,
            score: score,
            qstart: _qstart,
            qend: _qend,
            evalue: evaule,
            kick: false,
            frame: frame,
            sstart: sstart,
            send: send,
            seq: None,
            pident: pident,
            reference_hits: vec![ReferenceHit {
                target: ref_header,
                sstart,
                send,
            }],
            length: _qend as i32 - _qstart as i32 + 1,
            full_header: _full_header,
        }
    }


    fn convert_reference_hits(&self) -> Vec<String> {
        self.reference_hits
            .iter()
            .map(|refhit| refhit.to_json())
            .collect()
    }

    fn to_json(&self) -> String {
        let dict = serde_json::json!({
            "header": self.full_header,
            "seq": self.seq,
            "ref_taxon": self.reftaxon,
            "ali_start": self.sstart,
            "ali_end": self.qend,
            "reference_hits": self.convert_reference_hits()
        });
        serde_json::to_string(&dict).unwrap()
    }
}
