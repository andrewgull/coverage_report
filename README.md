Code to create coverage report with interactive plots.

## How to use

Edit 

- [input files](config/input_files.csv),
- [run IDs](data/run_IDs.tsv),

and add them to the params file (e.g. [params_v2.yml](params_v2.yml) or [params_v1.yml](params_v1.yml)) for the design version you want to use.

Then run:

```bash
pixi run render-v2
```

or

```bash
pixi run render-v1
```

### Note on parameters:

 - `gene_chunk_size` is the number of genes to show per page in the interactive plot.

 - `run_ids` is a tab-separated file with at least two columns: `run_ID` and `sample`, containing run IDs and sample names for each run (to find them check your samplesheet).
