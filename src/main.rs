use hansen::{Record, read_solvents};
fn main() {

    let solvs: Vec<Record> = read_solvents("solvents.csv".to_string());
    println!("{:?}", solvs.len())
}
