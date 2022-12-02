use hansen::{distance, line_segment, mixture, read_data, Drug, Solvent};
use std::thread;

fn main() {
    let drugs: Vec<Drug> = read_data::<Drug>("data/drug_list.csv".to_string());
    let mut handles = Vec::new();
    let solves: Vec<Solvent> = read_data::<Solvent>("data/solvents.csv".to_string());
    for drug in drugs {
        let solvs = solves.clone();
        let handle = thread::spawn(move || {
            println!("Starting Thread: {}", drug.drug);
            let mut closest: f32 = f32::MAX;
            let mut closest_a: String = "".to_string();
            let mut closest_b: String = "".to_string();
            for solvent_a in &solvs {
                let solvent_a: &Solvent = &solvent_a;
                for solvent_b in &solvs {
                    let solvent_b: &Solvent = &solvent_b;

                    let (start, end) = line_segment(solvent_a, solvent_b);

                    let c: f32 = distance(&drug, &start, &end);
                    if closest > c {
                        closest = c;
                    }
                }
            }

            println!("{}", closest);
        });
        handles.push(handle)
    }

    for handle in handles {
        handle.join().unwrap();
    }
}
