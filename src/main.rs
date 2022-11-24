use hansen::{distance, read_solvents, Record};

fn main() {
    let solvs: Vec<Record> = read_solvents("solvents.csv".to_string());
    let mut closest: f32 = f32::MAX;
    for solvent_a in &solvs {
        let solvent_a: &Record = solvent_a;
        for solvent_b in &solvs {
            let solvent_b: &Record = solvent_b;
            let c: f32 = distance(solvent_a, solvent_b);
            if closest < c {
                closest = c;
            }
        }
    }

    println!("{:?}", closest);
}
