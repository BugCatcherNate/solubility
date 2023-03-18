use clap::{arg, Command};
use hansen::{write_csv, TopN};

fn main() {
    let matches = Command::new("solvent_blend")
        .author("Nathan Thompson")
        .arg(
            arg!(--max_results <VALUE>)
                .help("Number of solvent blends to return")
                .required(true)
                .value_parser(clap::value_parser!(usize)),
        )
        .arg(
            arg!(--n <VALUE>)
                .help("Number of closest solvent blends to return for each drug")
                .required(true)
                .value_parser(clap::value_parser!(usize)),
        )
        .arg(
            arg!(--solvents <PATH>)
                .help("Path to csv list of solvents")
                .required(true)
                .value_parser(clap::value_parser!(String)),
        )
        .arg(
            arg!(--drugs <PATH>)
                .help("Path to csv list of drugs")
                .required(true)
                .value_parser(clap::value_parser!(String)),
        )
        .arg(
            arg!(--output <PATH>)
                .help("Path to output results")
                .required(false)
                .default_value("blend_ratios.csv")
                .value_parser(clap::value_parser!(String)),
        )
        .get_matches();

    let max_results: usize = *matches.get_one::<usize>("max_results").expect("required");

    let n: usize = *matches.get_one::<usize>("n").expect("required");
    let drugs_file_path = matches.get_one::<String>("drugs").expect("required");
    let sovlents_file_path = matches.get_one::<String>("solvents").expect("required");
    let blend_ratios_path = matches.get_one::<String>("output").unwrap();
    let tn: TopN = TopN::new(
        n.to_owned(),
        drugs_file_path.to_owned(),
        sovlents_file_path.to_owned(),
        max_results,
    );
    let res = tn.calculate();
    write_csv(res.clone(), blend_ratios_path.to_string());
}
