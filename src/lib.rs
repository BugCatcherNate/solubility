use csv::Reader;
use nalgebra::Vector3;
use rayon::prelude::{IntoParallelIterator, ParallelIterator};
use std::collections::HashMap;
use std::time::Instant;

use serde::de::DeserializeOwned;
use serde::{Deserialize, Serialize};
use std::cmp::Ordering::Equal;
#[derive(Debug, Deserialize, Clone, Serialize)]
pub struct Solvent {
    pub id: i32,
    pub solvent: String,
    pub d_d: f32,
    pub d_p: f32,
    pub d_h: f32,
}

#[derive(Debug, Clone, Deserialize)]
pub struct Drug {
    pub id: i32,
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

#[derive(Debug, Clone, Serialize, Copy)]
pub struct Solution {
    pub mix_id: i32,
    pub drug_id: i32,
    pub distance: f32,
}

#[derive(Debug, Clone, Serialize)]
pub struct FinalSolution {
    pub drug: String,
    pub mix_id: i32,
    pub solvent_a: String,
    pub solvent_a_ratio: f32,
    pub solvent_b: String,
    pub solvent_b_ratio: f32,
    pub hansen_distance: f32,
}
pub fn distance(drug: &Drug, start: &Vector3<f32>, end: &Vector3<f32>) -> f32 {
    let drug_params: Vector3<f32> = Vector3::new(drug.d_d, drug.d_p, drug.d_h);
    let num: f32 = Vector3::norm(&(end - start).cross(&(start - drug_params)));
    let dom: f32 = Vector3::norm(&(end - start));

    let res = num / dom;
    res
}

pub fn cantor(a: i32, b: i32) -> i32 {
    let sum: i32 = (a + b + 1) * (a + b);
    let triangle_sum: i32 = sum / 2;
    triangle_sum + a
}

pub fn inv_cantor(c: i32) -> (i32, i32) {
    let n = ((-1.0 + ((8 * c + 1) as f64).sqrt()) / 2.0).floor() as i32;
    let a = c - ((n + 1) * n) / 2;
    let b = n - a;
    (a, b)
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
        if dist_diff <= last_diff {
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
        sol_params: Vector3::new(mix_d_d, mix_d_p, mix_d_h),
    };
    new_blend
}

pub fn line_segment(a: &Solvent, b: &Solvent) -> (Vector3<f32>, Vector3<f32>) {
    let start: SolventMix = mixture(a, b, 0.9);
    let end: SolventMix = mixture(a, b, 0.1);

    (start.sol_params, end.sol_params)
}

pub fn write_data(solution: Vec<Solution>, path: String) {
    let mut wrt = csv::Writer::from_path(path).unwrap();

    for sol in solution {
        wrt.serialize(sol).unwrap();
    }

    wrt.flush().unwrap();
}

pub fn write_hash(solution: Vec<(&i32, &i32)>, path: String) {
    let mut wrt = csv::Writer::from_path(path).unwrap();

    for sol in solution {
        wrt.serialize(sol).unwrap();
    }

    wrt.flush().unwrap();
}

pub fn write_results(solution: Vec<FinalSolution>, path: String) {
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

pub struct TopN {
    n: usize,
    drugs_file: String,
    solves_file: String,
    max_results: usize,
}
impl TopN {
    pub fn new(n: usize, drugs_file: String, solves_file: String, max_results: usize) -> Self {
        Self {
            n,
            drugs_file,
            solves_file,
            max_results,
        }
    }

    pub fn calculate(&self) -> Vec<FinalSolution> {
        let drugs = read_data::<Drug>(self.drugs_file.to_string());
        let solves = read_data::<Solvent>(self.solves_file.to_string());

        let max_capacity = 100000;

        let par_iter = drugs.clone().into_par_iter().map(|drug| {
            let solvs = solves.clone();
            let mut temp_capacity: usize = max_capacity.clone();
            let mut top_mixes: Vec<Solution> = Vec::new();
            println!("Starting Thread: {}", drug.drug);
            let start = Instant::now();
            for solvent_a in &solvs {
                let solvent_a: &Solvent = &solvent_a;
                for solvent_b in &solvs {
                    if solvent_a.id < solvent_b.id {
                        let solvent_b: &Solvent = &solvent_b;

                        let (start, end) = line_segment(solvent_a, solvent_b);
                        let temp_solvent_a = solvent_a.clone();
                        let temp_solvent_b = solvent_b.clone();
                        let c: f32 = distance(&drug, &start, &end);
                        let temp_solution = Solution {
                            mix_id: cantor(temp_solvent_a.id, temp_solvent_b.id),
                            drug_id: drug.id,
                            distance: c,
                        };

                        top_mixes.push(temp_solution);

                        if top_mixes.len() > temp_capacity {
                            top_mixes.sort_by(|a, b| {
                                a.distance.partial_cmp(&b.distance).unwrap_or(Equal)
                            });
                            top_mixes = top_mixes.split_at(temp_capacity).0.to_vec();

                            temp_capacity *= 2;
                        }
                    }
                }
            }

            top_mixes.sort_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap_or(Equal));

            let duration = start.elapsed();
            println!("Finished Thread: {} in {:?} ", drug.drug, duration);
            top_mixes.split_at(self.n).0.to_vec()
        });

        let mut counts: HashMap<i32, i32> = HashMap::new();
        let res: Vec<Solution> = par_iter.flatten().collect();

        for r in &res {
            let new_count = match counts.get(&r.mix_id) {
                Some(count) => count + 1,
                None => 1,
            };
            counts.insert(r.mix_id, new_count);
        }

        let mut count_vec: Vec<_> = counts.iter().collect();

        count_vec.sort_by(|a, b| b.1.cmp(a.1));
        let final_counts = count_vec.split_at(self.max_results).0.to_vec();
        write_hash(final_counts.clone(), "mix_counts.csv".to_string());

        let mut final_results: Vec<FinalSolution> = Vec::new();
        let mut final_mixes: Vec<i32> = Vec::new();
        for x in final_counts {
            final_mixes.push(*x.0);
        }

        let temp_mixes: Vec<Solution> = res
            .clone()
            .into_iter()
            .filter(|s| final_mixes.contains(&s.mix_id))
            .collect();

        for mix in temp_mixes {
            let (solv_a_id, solv_b_id) = inv_cantor(mix.mix_id);
            let solv_a: Solvent = solves
                .clone()
                .into_iter()
                .find(|s| s.id == solv_a_id)
                .unwrap();

            let solv_b: Solvent = solves
                .clone()
                .into_iter()
                .find(|s| s.id == solv_b_id)
                .unwrap();

            let drug: Drug = drugs
                .clone()
                .into_iter()
                .find(|s| s.id == mix.drug_id)
                .unwrap();

            let (x_a, x_b): (f32, f32) = mix_solver(&solv_a, &solv_b, &drug, mix.distance);
            let temp_res: FinalSolution = FinalSolution {
                drug: drug.clone().drug,
                mix_id: mix.mix_id,
                solvent_a: solv_a.solvent,
                solvent_b: solv_b.solvent,
                hansen_distance: mix.distance,
                solvent_a_ratio: x_a * 100.0,
                solvent_b_ratio: x_b * 100.0,
            };
            final_results.push(temp_res);
        }

        final_results
    }
}


#[cfg(test)]
mod tests {
    use crate::{cantor, inv_cantor, standard_dist};

    #[test]
    fn test_cantor() {
        let test_a = 12222;
        let test_b = 7;
        let c = cantor(test_a, test_b);
        let (a, b) = inv_cantor(c);
        assert_eq!(a, test_a);
        assert_eq!(b, test_b);
    }


    #[test]
    fn test_standard_dist() {
        let point_a = (1.0, 1.0, 1.0);
        let point_b= (2.0, 2.0,2.0);
        let dist = standard_dist(point_a.0, point_a.1, point_a.2, point_b.0, point_b.1, point_b.2);
        assert_eq!(1.0, dist);
    }
}
