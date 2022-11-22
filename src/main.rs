
use csv::{Reader};
use serde::Deserialize;
#[derive(Debug, Deserialize)]
struct Record {
    id: i32,
    solvent: String,
    dD: f32,
    dP: f32,
    dH: f32

}
fn main() {

    let mut reader = Reader::from_path("solvents.csv").unwrap();
    let mut solvents: Vec<Record> = Vec::new();
    for result in reader.deserialize() {
        let record: Record = result.unwrap();
        solvents.push(record);
    }

    println!("{:?}", solvents.len())
}
