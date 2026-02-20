Code to create coverage report with interactive plots.

## How to use

Edit [input files](config/input_files.csv) and [run IDs](data/run_IDs.tsv), then run

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
