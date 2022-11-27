use hansen::{distance, mixture, read_data, Drug, Solvent};
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

                    let new_mix = mixture(solvent_a, solvent_b, 0.5);

                    let c: f32 = distance(&new_mix, &drug);
                    if closest > c {
                        closest = c;
                        closest_a = new_mix.solvent_a;
                        closest_b = new_mix.solvent_b;
                    }
                }
            }

            println!("{},{}", closest_a, closest_b);
        });
        handles.push(handle)
    }

    for handle in handles {
        handle.join().unwrap();
    }
}
