use csv::Reader;
use nalgebra::{Vector3};
use serde::de::DeserializeOwned;
use serde::{Deserialize, Serialize};
#[derive(Debug, Deserialize, Clone, Serialize)]
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


#[derive(Debug, Clone, Serialize)]
pub struct Solution {
    pub mix_id : i32,
    pub solvent_a: i32,
    pub solvent_b: i32,
    pub distance: f32
}
pub fn distance(drug: &Drug, start: &Vector3<f32>, end: &Vector3<f32> ) -> f32 {
    let drug_params: Vector3<f32> = Vector3::new(drug.d_d, drug.d_p, drug.d_h);
    let num: f32 = Vector3::norm(&(end - start).cross(&(start - drug_params)));  
    let dom: f32 = Vector3::norm(&(end - start));

    let res = num/dom;
    res
}

pub fn cantor(a: i32, b: i32) -> i32 {


    let sum: i32 = a + b;
    let triangle_sum: i32 = sum * (sum + 1)/2; 
    triangle_sum + b

}

pub fn standard_dist(a_x: f32, a_y: f32, a_z: f32, b_x: f32, b_y: f32, b_z: f32) -> f32 {

        ((b_x - a_x).powi(2) + (b_y - a_y).powi(2) + (b_z - a_z).powi(2)).sqrt()


}

pub fn mix_solver(a: &Solvent, b: &Solvent, drug: &Drug, dist: f32) -> (f32, f32) {
    let mut r_a: f32 = 0.9;
    let mut r_b: f32 = 1.0 - r_a;
    let mut last_diff = 1000000000.0;
    let mut best_r_a: f32 = 0.0;
    let mut best_r_b: f32 = 0.0;
    while r_a >= 0.1 {
        let a_x: f32 = r_a * a.d_d + r_b * b.d_d;
        let a_y: f32 = r_a * a.d_p + r_b * b.d_p;
        let a_z: f32 = r_a * a.d_h + r_b * b.d_h;
        let temp_dist = standard_dist(a_x, a_y, a_z, drug.d_d, drug.d_p, drug.d_h);
        let dist_diff = (temp_dist - dist).abs();
        if dist_diff <= last_diff{
            last_diff = dist_diff; 
            best_r_a = r_a;
            best_r_b = r_b;
        }
        r_a -= 0.01;
        r_b = 1.0 - r_a;
        }

    assert!(best_r_a != 0.0);

    assert!(best_r_b != 0.0);

    (best_r_a, best_r_b)

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

pub fn write_data(solution: Vec<Solution>, path: String){

    let mut wrt = csv::Writer::from_path(path).unwrap();

    for sol in solution {

        wrt.serialize(sol).unwrap();

    }

    wrt.flush().unwrap();

}

pub fn write_hash(solution: Vec<(&i32,&i32)>, path: String){

    let mut wrt = csv::Writer::from_path(path).unwrap();

    for sol in solution {

        wrt.serialize(sol).unwrap();

    }

    wrt.flush().unwrap();

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
