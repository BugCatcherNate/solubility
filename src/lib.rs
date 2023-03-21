use csv::Reader;

use log::{debug, info};
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
    // this function returns a f32 value representing the shortest distance between the Drug and the line segment defined by the start and end points.
    // this function takes in three arguments:
    // reference to a Drug struct
    // reference to a Vector3<f32> representing the start point
    // reference to another Vector3<f32> representing the end point.

    let res: f32;
    let drug_params: Vector3<f32> = Vector3::new(drug.d_d, drug.d_p, drug.d_h);
    let a_b = end - start;
    let b_e = drug_params - end;
    let a_e = drug_params - start;
    let ab_be = a_b.dot(&b_e);
    let ab_ae = a_b.dot(&a_e);

    if ab_be > 0.0 {
        // end point is closer to the Drug than the start point
        res = standard_dist(
            end.x,
            end.y,
            end.z,
            drug_params.x,
            drug_params.y,
            drug_params.z,
        );
    } else if ab_ae < 0.0 {
        // start point is closer to the Drug than the end point
        res = standard_dist(
            start.x,
            start.y,
            start.z,
            drug_params.x,
            drug_params.y,
            drug_params.z,
        );
    } else {
        // calculate perpendicular distance
        let num: f32 = Vector3::norm(&(end - start).cross(&(start - drug_params)));
        let dom: f32 = Vector3::norm(&(end - start));

        res = num / dom;
    };

    res
}

pub fn cantor(a: i32, b: i32) -> i32 {
    // this function returns the cantor pairing of integers a and b
    // cantor pairing is used to generate a unique identifier for a combination of two integers

    let sum: i32 = (a + b + 1) * (a + b);
    let triangle_sum: i32 = sum / 2;
    let res = triangle_sum + a;
    res
}

pub fn inv_cantor(c: i32) -> (i32, i32) {
    // this function reversed the cantor paring
    // result is the original integers used to create the unique identifier
    let n = ((-1.0 + ((8 * c + 1) as f64).sqrt()) / 2.0).floor() as i32;
    let a = c - ((n + 1) * n) / 2;
    let b = n - a;
    (a, b)
}

pub fn standard_dist(a_x: f32, a_y: f32, a_z: f32, b_x: f32, b_y: f32, b_z: f32) -> f32 {
    // standard Euclidean distance of 2 3d points
    ((b_x - a_x).powi(2) + (b_y - a_y).powi(2) + (b_z - a_z).powi(2)).sqrt()
}

pub fn mix_solver(a: &Solvent, b: &Solvent, drug: &Drug, dist: f32) -> (f32, f32) {
    let b_params = Vector3::new(b.d_d, b.d_p, b.d_h);
    let mut last_diff = f32::INFINITY;
    let mut best_r_a = 0.0;

    let mut best_dist: f32 = f32::INFINITY;
    debug!{"mix solver for {}, {} for distance {}", a.solvent, b.solvent, dist};
    for r_a in (10..=90).step_by(1).map(|r| r as f32 / 100.0) {
        let r_b = 1.0 - r_a;
        let a_params = Vector3::new(a.d_d, a.d_p, a.d_h);
        let mixed_params = a_params.scale(r_a) + b_params.scale(r_b);
        let temp_dist = standard_dist(
            mixed_params.x,
            mixed_params.y,
            mixed_params.z,
            drug.d_d,
            drug.d_p,
            drug.d_h,
        );
        let dist_diff = (temp_dist - dist).abs();

        if dist_diff <= last_diff {
            last_diff = dist_diff;
            best_r_a = r_a;
            best_dist = temp_dist;
        }
    }

    debug!{"mix solver for {}, {} for distance {}. Found ratio {}, {} with distnce {}", a.solvent, b.solvent, dist, best_r_a, 1.0 - best_r_a, best_dist};
    (best_r_a, 1.0 - best_r_a)
}

pub fn mixture(a: &Solvent, b: &Solvent, r_a: f32) -> SolventMix {
    let r_b = 1.0 - r_a;
    let sol_params = Vector3::new(
        r_a * a.d_d + r_b * b.d_d,
        r_a * a.d_p + r_b * b.d_p,
        r_a * a.d_h + r_b * b.d_h,
    );
    SolventMix {
        solvent_a: a.solvent.clone(),
        solvent_b: b.solvent.clone(),
        ratio_a: r_a,
        ratio_b: r_b,
        sol_params,
    }
}

//TODO add tests
pub fn line_segment(a: &Solvent, b: &Solvent) -> (Vector3<f32>, Vector3<f32>) {
    let start: SolventMix = mixture(a, b, 0.9);
    let end: SolventMix = mixture(a, b, 0.1);

    (start.sol_params, end.sol_params)
}

pub fn write_csv<T: serde::Serialize>(data: Vec<T>, path: String) {

    info!("writing to: {}", path);
    let mut wrt = csv::Writer::from_path(path).unwrap();

    for item in data {
        wrt.serialize(item).unwrap();
    }

    wrt.flush().unwrap();
}

pub fn read_csv<'a, T: DeserializeOwned>(path: String) -> Vec<T> {

    info!("reading from: {}", path);
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
        let mut drugs = read_csv::<Drug>(self.drugs_file.to_string());
        let mut solves = read_csv::<Solvent>(self.solves_file.to_string());
        let drugs_before_doubled = drugs[0].d_d;
        let solves_before_doubled = solves[0].d_d;
        for drug in &mut drugs  {
              drug.d_d *= 2.0;
        }
        for solvent in &mut solves  {
              solvent.d_d *= 2.0;
        }
        assert!((drugs[0].d_d == drugs_before_doubled * 2.0) & (solves[0].d_d == solves_before_doubled * 2.0)); 
        let max_capacity = 100000;

        let par_iter = drugs.clone().into_par_iter().map(|drug| {
            let solvs = solves.clone();
            let mut temp_capacity: usize = max_capacity.clone();
            let solvs_len = solvs.len();
            let mut top_mixes: Vec<Solution> = Vec::with_capacity(solvs_len * (solvs_len - 1) / 2);

            info!("Starting Thread: {}", drug.drug);
            let start = Instant::now();
            for (i, solvent_a) in solvs.iter().enumerate() {
                for solvent_b in solvs.iter().skip(i + 1) {
                    let (start, end) = line_segment(solvent_a, solvent_b);
                    let c: f32 = distance(&drug, &start, &end);
                    let temp_solution = Solution {
                        mix_id: cantor(solvent_a.id, solvent_b.id),
                        drug_id: drug.id,
                        distance: c,
                    };

                    top_mixes.push(temp_solution);

                    if top_mixes.len() > temp_capacity {
                        top_mixes
                            .sort_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap_or(Equal));
                        top_mixes = top_mixes.split_at(temp_capacity).0.to_vec();

                        temp_capacity *= 2;
                    }
                }
            }

            top_mixes.sort_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap_or(Equal));

            let duration = start.elapsed();

            info!("Finished Thread: {} in {:?} ", drug.drug, duration);
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
        write_csv(final_counts.clone(), "mix_counts.csv".to_string());

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

            let (x_a, x_b) = mix_solver(&solv_a, &solv_b, &drug, mix.distance);
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
    use float_cmp::approx_eq;
    use nalgebra::Vector3;

    // use crate::{cantor, distance, inv_cantor, standard_dist, Drug};
    use super::*;

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
        let point_b = (2.0, 2.0, 2.0);
        let dist = standard_dist(
            point_a.0, point_a.1, point_a.2, point_b.0, point_b.1, point_b.2,
        );
        assert!(approx_eq!(f32, 1.732, dist, epsilon = 0.001));
    }

    #[test]
    fn test_dist_a() {
        let start = Vector3::new(2.0, 2.0, 2.0);
        let end = Vector3::new(3.0, 3.0, 3.0);
        let test_drug = Drug {
            id: 1,
            drug: "test_drug".to_string(),
            d_d: 1.0,
            d_p: 1.0,
            d_h: 1.0,
        };
        let d = distance(&test_drug, &start, &end);

        assert!(approx_eq!(f32, 1.73205, d, epsilon = 0.001));
    }
    #[test]
    fn test_dist_b() {
        let start = Vector3::new(2.0, 2.0, 2.0);
        let end = Vector3::new(3.0, 3.0, 3.0);
        let test_drug = Drug {
            id: 1,
            drug: "test_drug".to_string(),
            d_d: 5.0,
            d_p: 5.0,
            d_h: 5.0,
        };
        let d = distance(&test_drug, &start, &end);

        assert!(approx_eq!(f32, 3.464102, d, epsilon = 0.001));
    }

    #[test]
    fn test_dist_c() {
        let start = Vector3::new(5.0, 2.0, 1.0);
        let end = Vector3::new(3.0, 1.0, -1.0);
        let test_drug = Drug {
            id: 1,
            drug: "test_drug".to_string(),
            d_d: 0.0,
            d_p: 2.0,
            d_h: 3.0,
        };
        let d = distance(&test_drug, &start, &end);
        println!("{}", d);

        assert!(approx_eq!(f32, 5.0, d, epsilon = 0.001));
    }

    #[test]
    fn test_mixture() {
        let a = Solvent {
            id: 1,
            solvent: "A".to_string(),
            d_d: 1.0,
            d_p: 2.0,
            d_h: 3.0,
        };
        let b = Solvent {
            id: 2,
            solvent: "B".to_string(),
            d_d: 4.0,
            d_p: 5.0,
            d_h: 6.0,
        };
        let r_a = 0.6;
        let expected = SolventMix {
            solvent_a: "A".to_string(),
            solvent_b: "B".to_string(),
            ratio_a: 0.6,
            ratio_b: 0.4,
            sol_params: Vector3::new(2.2, 3.2, 4.2),
        };
        let result = mixture(&a, &b, r_a);
        assert_eq!(result.solvent_a, expected.solvent_a);
        assert_eq!(result.solvent_b, expected.solvent_b);
        assert!(approx_eq!(
            f32,
            result.ratio_a,
            expected.ratio_a,
            epsilon = 0.001
        ));
        assert!(approx_eq!(
            f32,
            result.ratio_b,
            expected.ratio_b,
            epsilon = 0.001
        ));
        assert!(approx_eq!(
            f32,
            result.sol_params.x,
            expected.sol_params.x,
            epsilon = 0.001
        ));
        assert!(approx_eq!(
            f32,
            result.sol_params.y,
            expected.sol_params.y,
            epsilon = 0.001
        ));
        assert!(approx_eq!(
            f32,
            result.sol_params.z,
            expected.sol_params.z,
            epsilon = 0.001
        ));
    }

    #[test]
    fn test_mix_solver_a() {
        let a = Solvent {
            id: 1,
            solvent: "solvent1".to_string(),
            d_d: 1.0,
            d_p: 2.0,
            d_h: 3.0,
        };
        let b = Solvent {
            id: 2,
            solvent: "solvent2".to_string(),
            d_d: 4.0,
            d_p: 5.0,
            d_h: 6.0,
        };
        let drug = Drug {
            id: 1,
            drug: "drug".to_string(),
            d_d: 7.0,
            d_p: 8.0,
            d_h: 9.0,
        };
        let dist = 5.7157;

        let (r_a, r_b) = mix_solver(&a, &b, &drug, dist);
        assert_eq!(r_a, 0.1, "Expected r_a to be 0.1");
        assert_eq!(r_b, 0.9, "Expected r_b to be 0.9");
    }
    #[test]
    fn test_mix_solver_b() {
        let a = Solvent {
            id: 1,
            solvent: "solvent1".to_string(),
            d_d: 1.0,
            d_p: 2.0,
            d_h: 3.0,
        };
        let b = Solvent {
            id: 2,
            solvent: "solvent2".to_string(),
            d_d: 4.0,
            d_p: 5.0,
            d_h: 6.0,
        };
        let drug = Drug {
            id: 1,
            drug: "drug".to_string(),
            d_d: 7.0,
            d_p: 8.0,
            d_h: 9.0,
        };
        let dist = 6.02;

        let (r_a, r_b) = mix_solver(&a, &b, &drug, dist);

        assert!(approx_eq!(f32, r_a, 0.16, epsilon = 0.001));
        assert!(approx_eq!(f32, r_b, 0.84, epsilon = 0.001));
    }
}
