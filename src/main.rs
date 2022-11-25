use hansen::{distance, read_solvents, Record, mixture};

fn main() {
    let solvs: Vec<Record> = read_solvents("solvents.csv".to_string());
    let mut closest: f32 = f32::MAX;
    let test_drug = Record(id:0, solvent:"meme", d_d: 1.0, d_p: 1.0, d_h: 1.0);
    for solvent_a in &solvs {
        let solvent_a: &Record = solvent_a;
        for solvent_b in &solvs {
            let solvent_b: &Record = solvent_b;

            let new_mix = mixture(solvent_a, solvent_b;

            let c: f32 = distance(new_mix, solvent_b);
            if closest < c {
                closest = c;
            }
        }
    }

    println!("{:?}", closest);
}
