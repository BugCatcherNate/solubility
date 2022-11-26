use csv::Reader;
use serde::Deserialize;
#[derive(Debug, Deserialize)]
pub struct Record {
    pub id: i32,
    pub solvent: String,
    pub d_d: f32,
    pub d_p: f32,
    pub d_h: f32,
}

pub struct Solvent_mix<'solvent> {
    pub solvent_a: &'solvent Record,
    pub solvent_b: &'solvent Record,
    pub ratio_a: f32,
    pub ratio_b: f32,
    pub d_d: f32,
    pub d_p: f32,
    pub d_h: f32,
}

pub fn distance(a: &Solvent_mix, b: &Record) -> f32 {
    let num = 4.0 * (b.d_d - a.d_d).powi(2) + (b.d_p - a.d_p).powi(2) + (b.d_h - a.d_h).powi(2);
    let res = num.sqrt();
    res
}

pub fn mixture<'solvent>(a: &'solvent Record, b: &'solvent Record, r_a: f32) -> Solvent_mix<'solvent>{
    let r_b: f32 = 1.0 - r_a;
    let mix_d_d :f32 = r_a * a.d_d + r_b * b.d_d; 
    let mix_d_p :f32 = r_a * a.d_p + r_b * b.d_p; 
    let mix_d_h :f32 = r_a * a.d_h + r_b * b.d_h; 
    let new_blend = Solvent_mix{solvent_a: a, solvent_b: b, ratio_a: r_a, ratio_b: r_b, d_d: mix_d_d, d_p: mix_d_p, d_h: mix_d_h}; 
    new_blend
}



pub fn read_solvents(path: String) -> Vec<Record> {
    let mut reader = Reader::from_path(path).unwrap();
    let mut solvents: Vec<Record> = Vec::new();
    for result in reader.deserialize() {
        let record: Record = result.unwrap();
        solvents.push(record);
    }
    solvents
}
