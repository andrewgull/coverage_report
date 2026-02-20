Code to create coverage report with interactive plots.

## How to use

Edit `params` section to set paths to [input files](config/input_files.csv) and [run IDs](data/run_IDs.tsv):

```yaml
params:
  input_files: "config/input_files.csv"
  gene_chunk_size: 20
  run_ids: "data/run_IDs.tsv"
```

`gene_chunk_size` is the number of genes to show per page in the interactive plot.

`run_ids` is a tab-separated file with at least two columns: `run_ID` and `sample`, containing run IDs and sample names for each run (to find them check your samplesheet).

To render the report use `pixi run render-v2` for the new version or `pixi run render-v1` for the old version.
