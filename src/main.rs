use hansen::{
    cantor, distance, line_segment, mix_solver, read_data, write_hash, write_results, Drug,
    FinalSolution, Solution, Solvent,
};
use nalgebra::base;
use rayon::prelude::*;
use std::cmp::Ordering::Equal;
use std::collections::HashMap;
use std::env;
use std::time::Instant;

fn main() {
    let mut counts: HashMap<i32, i32> = HashMap::new();
    let args: Vec<String> = env::args().collect();

    let max_results: usize = args[3].parse::<usize>().unwrap();
    let base_sol_a: i32 = args[1].parse::<i32>().unwrap();
    let base_sol_b: i32 = args[2].parse::<i32>().unwrap();
    let max_capacity: usize = 100000;

    let drugs: Vec<Drug> = read_data::<Drug>("data/drug_list.csv".to_string());
    let solves: Vec<Solvent> = read_data::<Solvent>("data/solvents.csv".to_string());


    let base_solv_a: Solvent = solves
                .clone()
                .into_iter()
                .find(|s| s.id == base_sol_a)
                .unwrap();

    let base_solv_b: Solvent = solves
                .clone()
                .into_iter()
                .find(|s| s.id == base_sol_b)
                .unwrap();

    println!("Comparing Solvent Mixes against {} and {}", base_solv_a.solvent, base_solv_b.solvent);
    let par_iter = drugs.into_par_iter().map(|drug| {

        let (base_start, base_end) = line_segment(&base_solv_a, &base_solv_b);
        
        let base_distance: f32 = distance(&drug, &base_start, &base_end);

        println!("Comparisoon Solvent Mixes {} and {} distance to {} is {:}", base_solv_a.solvent, base_solv_b.solvent, drug.drug, base_distance);
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
                        solvent_a: temp_solvent_a.id,
                        solvent_b: temp_solvent_b.id,
                        distance: c,
                    };

                    if c <= base_distance {

                        top_mixes.push(temp_solution);
                    }


//                    if top_mixes.is_empty() || top_mixes.len() < temp_capacity {
                        //top_mixes.push(temp_solution);
                        //if top_mixes.len() == temp_capacity {
                            //top_mixes.sort_by(|a, b| {
                                //a.distance.partial_cmp(&b.distance).unwrap_or(Equal)
                            //});
                        //}
                    //} else if top_mixes.last().unwrap().distance > c {
                        //top_mixes.push(temp_solution);

                        //top_mixes.par_sort_unstable_by(|a, b| {
                            //a.distance.partial_cmp(&b.distance).unwrap_or(Equal)
                        //});

                        //if top_mixes.len() > temp_capacity {
                            //top_mixes.pop();
                        //}

                        //temp_capacity *= 2;
                    //}
                }
            }
        }

        let final_mixes = top_mixes.split_at(max_results).0.to_vec();
        let mut final_results: Vec<FinalSolution> = Vec::new();
        for mix in &final_mixes {
            let solv_a: Solvent = solvs
                .clone()
                .into_iter()
                .find(|s| s.id == mix.solvent_a)
                .unwrap();

            let solv_b: Solvent = solvs
                .clone()
                .into_iter()
                .find(|s| s.id == mix.solvent_b)
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
        let duration = start.elapsed();
        println!("Finished Thread: {} in {:?} ", drug.drug, duration);
        final_results
    });
    let res: Vec<FinalSolution> = par_iter.flatten().collect();

    write_results(res.clone(), "results.csv".to_string());

    for r in res {
        let new_count = match counts.get(&r.mix_id) {
            Some(count) => count + 1,
            None => 1,
        };
        counts.insert(r.mix_id, new_count);
    }
    let mut count_vec: Vec<_> = counts.iter().collect();
    count_vec.sort_by(|a, b| b.1.cmp(a.1));
    write_hash(count_vec.split_at(50).0.to_vec(), "res.csv".to_string());
}
