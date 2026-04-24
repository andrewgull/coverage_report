# Chromosome Coverage Comparison Report — Design Spec

**Date:** 2026-04-24
**Status:** Approved

## Goal

A new standalone Quarto HTML report that renders per-chromosome coverage boxplots, showing a side-by-side comparison between two probe panel designs (v1 and v2). Each chromosome gets a single combined Plotly figure with one box per design, sharing a common y-axis scale.

---

## Parameters

File: `params_comparison.yml`

```yaml
bed_dir_v1: "/path/to/design1/mosdepth_output"
bed_dir_v2: "/path/to/design2/mosdepth_output"
```

Both parameters are paths to directories containing per-sample mosdepth `.regions.bed.gz` files — the same structure as the existing `mosdepth_dir` parameter in `params_v1.yml` / `params_v2.yml`.

---

## Data Loading

In the `setup` chunk of `chr_comparison_report.qmd`:

1. Call `list_bed_files(params$bed_dir_v1, ".regions.bed.gz")` and `list_bed_files(params$bed_dir_v2, ".regions.bed.gz")` to discover files in each directory.
2. Map `read_bed()` over each list to produce two long-format tibbles.
3. Add a `design` column (`"v1"` / `"v2"`) to each tibble.
4. `bind_rows()` the two tibbles into a single data frame with columns: `X1` (chr), `X2` (start), `X3` (end), `X4` (region), `X5` (coverage depth), `sample`, `design`.
5. Derive `chromosomes` dynamically from `unique(bed$X1)`, sorted with base `sort()` (produces chr1, chr10 … chr2 … chrX, chrY — alphabetical order is acceptable; no extra dependency needed).

No GTF annotation or bedtools shell chunk is needed — the chromosome name comes directly from `X1`.

---

## New Function: `make_chr_comparison_plot`

Location: `R/functions.R`

```
make_chr_comparison_plot(chromosome, bed)
```

- **Inputs:** `chromosome` (single character string, e.g. `"chr1"`), `bed` (combined data frame from both designs).
- **Behaviour:**
  - Filters `bed` to `X1 == chromosome`.
  - Warns and returns an empty figure if no rows match.
  - Calls `plot_ly()` with `x = ~design`, `y = ~X5`, `color = ~design`, `type = "box"`, `boxpoints = "all"`, `jitter = 0.3`, `pointpos = 0` — consistent with existing plot style in `functions.R`.
  - Layout: title `"<chromosome> depth of coverage"`, y-axis label `"depth"`, x-axis label `"design"`.
- **Returns:** A plotly figure object.

---

## Report Structure: `chr_comparison_report.qmd`

### YAML front matter

```yaml
title: "Chromosome Coverage Comparison"
subtitle: "Design v1 vs v2"
author: "A.G."
format:
  html:
    self-contained: true
    theme: "cosmo"
    toc: true
    toc-location: left
params:
  bed_dir_v1: "/path/to/v1"
  bed_dir_v2: "/path/to/v2"
```

### Chunks

| Chunk | Purpose |
|---|---|
| `setup` | `source(here("R", "functions.R"))`, read + bind both dirs, derive `chromosomes` vector |
| `plots` | `for` loop over `chromosomes`: print `##` heading + `make_chr_comparison_plot()` |

No bedtools shell chunk — no GTF/probe annotation required.

---

## Build Task

Add to `pixi.toml`:

```toml
render-comparison = "quarto render chr_comparison_report.qmd --output \"coverage_comparison_$(date +%Y-%m-%d_%H-%M-%S).html\" --execute-params params_comparison.yml"
```

---

## Files Changed / Created

| File | Action |
|---|---|
| `chr_comparison_report.qmd` | Create |
| `params_comparison.yml` | Create |
| `R/functions.R` | Add `make_chr_comparison_plot()` |
| `pixi.toml` | Add `render-comparison` task |

---

## Out of Scope

- Per-gene or per-sample breakdowns (covered by `main_report.qmd`)
- GTF annotation or bedtools intersection
- Filtering to canonical chromosomes only (all chromosomes in the data are shown)
