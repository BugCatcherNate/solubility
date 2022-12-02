use hansen::{distance, line_segment, read_data, Solution, Drug, Solvent};
use std::env;
use std::cmp::Ordering::Equal;
use std::thread;

fn main() {
    let args: Vec<String> = env::args().collect();
    let max_results: usize = args[1].parse::<usize>().unwrap();
    let drugs: Vec<Drug> = read_data::<Drug>("data/drug_list.csv".to_string());
    let mut handles = Vec::new();
    let solves: Vec<Solvent> = read_data::<Solvent>("data/solvents.csv".to_string());
    for drug in drugs {

        let solvs = solves.clone();
        let handle = thread::spawn(move || {

            let mut top_mixes: Vec<Solution> = Vec::with_capacity(max_results);
            println!("Starting Thread: {}", drug.drug);
            for solvent_a in &solvs {
                let solvent_a: &Solvent = &solvent_a;
                for solvent_b in &solvs {
                    let solvent_b: &Solvent = &solvent_b;

                    let (start, end) = line_segment(solvent_a, solvent_b);
                    let temp_solvent_a = solvent_a.clone();
                    let temp_solvent_b = solvent_b.clone();
                    let c: f32 = distance(&drug, &start, &end);
                    let temp_solution = Solution{ solvent_a: temp_solvent_a.solvent, solvent_b: temp_solvent_b.solvent, distance: c};
                    if top_mixes.is_empty(){
                    top_mixes.push(temp_solution);
                    } else if top_mixes.last().unwrap().distance > c {
                        top_mixes.push(temp_solution);
                        top_mixes.sort_by(|a,b| a.distance.partial_cmp(&b.distance).unwrap_or(Equal));
                       
                        if top_mixes.len() > max_results {
                            top_mixes.pop();
                        }
                    }

            

                }
            }

            println!("{:?}", top_mixes);
        });
        handles.push(handle)
    }

    for handle in handles {
        handle.join().unwrap();
    }
}
