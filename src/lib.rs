
use csv::{Reader};
use serde::Deserialize;
#[derive(Debug, Deserialize)]
pub struct Record {
   pub id: i32,
    pub solvent: String,
   pub  d_d: f32,
    pub d_p: f32,
   pub  d_h: f32

}


pub fn distance(a: &Record, b: &Record) -> f32{

        
    let num = (b.d_d-a.d_d).powi(2) + (b.d_p-a.d_p).powi(2) + (b.d_h-a.d_h).powi(2);  
    let res = num.sqrt();
    res
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

