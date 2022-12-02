use csv::Reader;
use nalgebra::{Vector3};
use serde::de::DeserializeOwned;
use serde::Deserialize;
#[derive(Debug, Deserialize, Clone)]
pub struct Solvent {
    pub id: i32,
    pub solvent: String,
    pub d_d: f32,
    pub d_p: f32,
    pub d_h: f32,
}

#[derive(Debug, Deserialize)]
pub struct Drug {
    pub drug: String,
    pub d_d: f32,
    pub d_p: f32,
    pub d_h: f32,
}

#[derive(Debug, Clone)]
pub struct SolventMix {
    pub solvent_a: String,
    pub solvent_b: String,
    pub ratio_a: f32,
    pub ratio_b: f32,
    pub sol_params: Vector3<f32>,
}

pub fn distance(drug: &Drug, start: &Vector3<f32>, end: &Vector3<f32> ) -> f32 {
    let drug_params: Vector3<f32> = Vector3::new(drug.d_d, drug.d_p, drug.d_h);
    let num: f32 = Vector3::norm(&(end - start).cross(&(start - drug_params)));  
    let dom: f32 = Vector3::norm(&(end - start));

    let res = num/dom;
    res
}

pub fn mixture(a: &Solvent, b: &Solvent, r_a: f32) -> SolventMix {
    let r_b: f32 = 1.0 - r_a;
    let mix_d_d: f32 = r_a * a.d_d + r_b * b.d_d;
    let mix_d_p: f32 = r_a * a.d_p + r_b * b.d_p;
    let mix_d_h: f32 = r_a * a.d_h + r_b * b.d_h;
    let new_blend = SolventMix {
        solvent_a: a.solvent.clone(),
        solvent_b: b.solvent.clone(),
        ratio_a: r_a,
        ratio_b: r_b,
        sol_params: Vector3::new(mix_d_d, mix_d_p, mix_d_h)
    };
    new_blend
}

pub fn line_segment(a: &Solvent, b: &Solvent) -> (Vector3<f32>, Vector3<f32>) {


    let start: SolventMix = mixture(a, b, 0.9);
    let end: SolventMix = mixture(a, b, 0.1);
    
    (start.sol_params, end.sol_params)

}

pub fn read_data<'a, T: DeserializeOwned>(path: String) -> Vec<T> {
    let mut reader = Reader::from_path(path).unwrap();
    let mut drugs: Vec<T> = Vec::new();
    for result in reader.deserialize() {
        let record: T = result.unwrap();
        drugs.push(record);
    }
    drugs
}
