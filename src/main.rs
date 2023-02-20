use hansen::{write_results, TopN};
use std::env;

//TODO add user help function
//TODO fancier input with tags
//TODO error handling
fn main() {
    let args: Vec<String> = env::args().collect();

    let max_results: usize = args[1].parse::<usize>().unwrap();

            let n: usize = args[2].parse::<usize>().unwrap();

            let tn: TopN = TopN::new(
                n.to_owned(),
                "data/drug_list.csv".to_owned(),
                "data/solvents.csv".to_owned(),
                max_results
            );
            let res = tn.calculate();
            write_results(res.clone(), "mix_ratios.csv".to_string());


    }

    
