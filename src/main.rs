use hansen::{distance, mixture, read_drugs, read_solvents, Drug, Record};
use std::thread;

fn main() {
    let drugs: Vec<Drug> = read_drugs("drug_list.csv".to_string());
    let mut handles = Vec::new();
    for drug in drugs {
        let handle = thread::spawn(move || {
            let solvs: Vec<Record> = read_solvents("solvents.csv".to_string());
            println!("Starting Thread: {}", drug.drug);
            let mut closest: f32 = f32::MAX;
            for solvent_a in &solvs {
                let solvent_a: &Record = &solvent_a;
                for solvent_b in &solvs {
                    let solvent_b: &Record = &solvent_b;

                    let new_mix = mixture(solvent_a, solvent_b, 0.5);

                    let c: f32 = distance(&new_mix, &drug);
                    if closest > c {
                        closest = c;
                    }
                }
            }

            println!("{:}", closest);
        });
        handles.push(handle)
    }

    for handle in handles {
        handle.join().unwrap();
    }
}
