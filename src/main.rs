use hansen::{write_results, BetterSolvent, TopN};
use std::env;

fn main() {
    let args: Vec<String> = env::args().collect();
    let mode: char = args[1].parse::<char>().unwrap();

    let max_results: usize = args[2].parse::<usize>().unwrap();
    match mode {
        's' => {
            let base_sol_a: i32 = args[3].parse::<i32>().unwrap();
            let base_sol_b: i32 = args[4].parse::<i32>().unwrap();

            let bs = BetterSolvent::new(
                base_sol_a,
                base_sol_b,
                "data/drug_list.csv".to_owned(),
                "data/solvents.csv".to_owned(),
            );
            let res = bs.calculate();
        
        }

        'n' => {
            let n: usize = args[2].parse::<usize>().unwrap();

            let tn: TopN = TopN::new(
                n.to_owned(),
                "data/drug_list.csv".to_owned(),
                "data/solvents.csv".to_owned(),
                max_results
            );
            let res = tn.calculate();
            write_results(res.clone(), "results.csv".to_string());

        }

        _ => {
            panic!("Option not recognized")
        }
    }

    
    }
