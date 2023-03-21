# Hansen Solubility

This is a rust application created for finding the top combinations of two component solvent blends for dissolving a list of specified drugs.
## Build From Source
```bash
cd solubility
cargo build --release
```
## Run Tests With

```bash
cd solubility
cargo test
```

## Usage

```bash
./target/release/hansen --solvents solvents.csv --drugs drugs.csv --n 100000 --max_results 100
```
***Input Files***

- solvents.csv

|id|solvent|d_d|d_p|d_h|
|--|-------|---|---|---|
|1|"1,1,2,2-Tetrachloroethane"|18.8|5.1|5.3|
|2|"1,1,2-Trichloroethane"|18.2|5.3|6.8|

- drugs.csv

|id|drug|d_d|d_p|d_h|
|--|-------|---|---|---|
|1|6-monoacetylmorphine|19.27|3.01|7.6|
|2|7-Aminoclonazepam|21.04|13.16|7.11|
|3|Alpha OH Alprazolam|20.68|9.73|5.96|

***Parameters***

      --solvents <PATH>      Path to csv list of solvents
      --drugs <PATH>         Path to csv list of drugs
      --max_results <VALUE>  Number of solvent blends to return [default: 100]
      --n <VALUE>            Number of closest solvent blends to return for each drug [default: 100000]
      --output <PATH>        Path to output results [default: blend_ratios.csv]
      -h, --help                 Print help

The output will be 2 csv files
- blend_ratios.csv
- mix_counts.csv

The blend_ratios.csv will consist of the target drug, the 2 solvents and their ratios needed for the blend and their distance to the drug in hansen space. The ID column is the unique ID of the solvent blend. It can be used to match it to the other output file (mix_counts.csv)

|drug|mix_id|solvent_a|solvent_a_ratio|solvent_b|solvent_b_ratio|hansen_distance|
|----|----|----|----|------|------|---------|
|6-monoacetylmorphine|19805600|1-Hydroxy-2-Naphthoic Acid|45.0|a-Methyl Styrene|55.0|0.58334|
|6-monoacetylmorphine|31457513|"2-Naphthalenecarboxylic Acid, 3-Hydroxy-"|45.0|a-Methyl Styrene|55.0|0.58334|

The mix_counts.csv will consist of the solvent blend ID and the number of solvents that it was in the top blends for based on the parameter specified by --n

|Blend ID| count|
|-------|-------|
|2935058|9|
|26358408|9|
## License

[MIT](https://choosealicense.com/licenses/mit/)
