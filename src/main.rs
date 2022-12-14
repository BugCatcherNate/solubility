use hansen::{
    cantor, distance, line_segment, read_data, write_data, write_hash, Drug, Solution, Solvent,
};
use rayon::prelude::*;
use std::cmp::Ordering::Equal;
use std::collections::HashMap;
use std::env;
use std::time::{Duration, Instant};

fn main() {
    let mut counts: HashMap<i32, i32> = HashMap::new();
    let args: Vec<String> = env::args().collect();
    let max_results: usize = args[1].parse::<usize>().unwrap();
    let drugs: Vec<Drug> = read_data::<Drug>("data/drug_list.csv".to_string());
    let solves: Vec<Solvent> = read_data::<Solvent>("data/solvents.csv".to_string());
    let par_iter = drugs.into_par_iter().map(|drug| {
        let solvs = solves.clone();

        let mut top_mixes: Vec<Solution> = Vec::with_capacity(10000);
        println!("Starting Thread: {}", drug.drug);
        let start = Instant::now();
        for solvent_a in &solvs {
            println!("{}, {}", drug.drug, solvent_a.id.to_string());
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
                        solvent_a: temp_solvent_a.solvent,
                        solvent_b: temp_solvent_b.solvent,
                        distance: c,
                    };
                    if top_mixes.is_empty() || top_mixes.len() < max_results {
                        top_mixes.push(temp_solution);
                        if top_mixes.len() == max_results {
                            top_mixes.sort_by(|a, b| {
                                a.distance.partial_cmp(&b.distance).unwrap_or(Equal)
                            });
                        }
                    } else if top_mixes.last().unwrap().distance > c {
                        top_mixes.push(temp_solution);

      top_mixes.par_sort_unstable_by(|a, b| {
                                a.distance.partial_cmp(&b.distance).unwrap_or(Equal)
                            });
 
                        if top_mixes.len() > max_results {
                            top_mixes.pop();
                        }
                    }
                }
            }
        }

        let duration = start.elapsed();
        println!("Finished Thread: {} in {:?} ", drug.drug, duration);
        top_mixes.split_at(max_results).0.to_vec()
    });

    let res: Vec<Solution> = par_iter.flatten().collect();

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
