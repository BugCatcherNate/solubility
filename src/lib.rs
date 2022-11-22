
use csv::{Reader};
use serde::Deserialize;
#[derive(Deserialize)]
pub struct Record {
    id: i32,
    solvent: String,
    d_d: f32,
    d_p: f32,
    d_h: f32

}

pub fn read_solvents(path: String) -> Vec<Record>{

    let mut reader = Reader::from_path(path).unwrap();
    let mut solvents: Vec<Record> = Vec::new();
    for result in reader.deserialize() {
        let record: Record = result.unwrap();
        solvents.push(record);
    }
    solvents
}
