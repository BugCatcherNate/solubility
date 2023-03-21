# Hansen Solubility

This is a rust application created for finding the top combinations of two component solvent blends for dissolving a list of specified drugs.
## Build From Source
```bash
cd solubility
cargo build --release
```
## Usage

```bash
./target/release/hansen --solvents solvents.csv --drugs drugs.csv --n 100000 --max_results 100
```
***Parameters***

      --solvents <PATH>      Path to csv list of solvents
      --drugs <PATH>         Path to csv list of drugs
      --max_results <VALUE>  Number of solvent blends to return [default: 100]
      --n <VALUE>            Number of closest solvent blends to return for each drug [default: 100000]
      --output <PATH>        Path to output results [default: blend_ratios.csv]
      -h, --help                 Print help

## License

[MIT](https://choosealicense.com/licenses/mit/)
