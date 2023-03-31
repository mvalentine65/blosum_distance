use pyo3::prelude::*;
use serde::Deserialize;
use serde::Serialize;
use serde_json::to_string;

#[pyclass]
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct ReferenceHit {
    #[pyo3(get, set)]
    target: String,
    #[pyo3(get, set)]
    sstart: u16,
    #[pyo3(get, set)]
    send: u16,
}

#[pymethods]
impl ReferenceHit {
    #[new]
    fn new(target: String, sstart: u16, send: u16) -> Self {
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
    gene: String,
    #[pyo3(get, set)]
    frame: i8,
    #[pyo3(get, set)]
    evalue: f64,
    #[pyo3(get, set)]
    score: f32,
    #[pyo3(get, set)]
    qstart: u16,
    #[pyo3(get, set)]
    qend: u16,
    #[pyo3(get, set)]
    sstart: u16,
    #[pyo3(get, set)]
    send: u16,
    #[pyo3(get, set)]
    pident: f32,
    #[pyo3(get, set)]
    reftaxon: String,
    #[pyo3(get, set)]
    kick: bool,
    #[pyo3(get, set)]
    seq: Option<String>,
    #[pyo3(get, set)]
    length: u16,
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
        frame: String,
        evaule: String,
        score: String,
        qstart: String,
        qend: String,
        sstart: String,
        send: String,
        pident: String,
        gene: String,
        reftaxon: String,
    ) -> Self {
        let _qstart;
        let _qend;
        let _full_header;
        let frame = frame.parse::<i8>().unwrap();
        let sstart = sstart.parse::<u16>().unwrap();
        let send= send.parse::<u16>().unwrap();
        if frame >= 0 {
            // normal values
            _qstart = qstart.parse::<u16>().unwrap();
            _qend = qend.parse::<u16>().unwrap();
            _full_header = format!("{}|[translate({})]", header, frame);
        } else {
            // frame < 0, so reverse and mark as revcomp
            _qstart = qend.parse::<u16>().unwrap();
            _qend = qstart.parse::<u16>().unwrap();
            _full_header = format!("{}|[revcomp]:[translate({})]", header, frame.abs());
        };
        Self {
            header: header,
            target: ref_header.clone(),
            gene: gene,
            reftaxon: reftaxon,
            score: score.parse::<f32>().unwrap(),
            qstart: _qstart,
            qend: _qend,
            evalue: evaule.parse::<f64>().unwrap(),
            kick: false,
            frame: frame,
            sstart: sstart,
            send: send,
            seq: None,
            pident: pident.parse::<f32>().unwrap(),
            reference_hits: vec![ReferenceHit {
                target: ref_header,
                sstart,
                send,
            }],
            length: _qend - _qstart + 1,
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
