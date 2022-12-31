use hansen::{write_hash, write_results, BetterSolvent};
use std::collections::HashMap;
use std::env;

fn main() {
    let mut counts: HashMap<i32, i32> = HashMap::new();
    let args: Vec<String> = env::args().collect();
    let mode: char = args[1].parse::<char>().unwrap();
    match mode {

        'n' => {
            let base_sol_a: i32 = args[2].parse::<i32>().unwrap();
    let base_sol_b: i32 = args[3].parse::<i32>().unwrap();

    let base_sol_b: i32 = args[4].parse::<i32>().unwrap();
        }
    }

    let max_results: usize = args[3].parse::<usize>().unwrap();

    let bs = BetterSolvent::new(
        base_sol_a,
        base_sol_b,
        max_results,
        "data/drug_list.csv".to_owned(),
        "data/solvents.csv".to_owned(),
    );

    let res = bs.calculate();
    write_results(res.clone(), "results.csv".to_string());

    for r in res {
        let new_count = match counts.get(&r.mix_id) {
            Some(count) => count + 1,
            None => 1,
        };
        counts.insert(r.mix_id, new_count);
    }
    let mut count_vec: Vec<_> = counts.iter().collect();
    count_vec.sort_by(|a, b| b.1.cmp(a.1));
    write_hash(count_vec.split_at(50).0.to_vec(), "res.csv".to_string());
}
