# Chromosome Comparison Report Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a standalone Quarto HTML report that renders one combined Plotly boxplot per chromosome, comparing per-region coverage depth between two mosdepth experiment designs side by side.

**Architecture:** A new `chr_comparison_report.qmd` reads two mosdepth output directories (one per design), binds the data with a `design` label column, discovers chromosomes dynamically from `X1`, and renders one `make_chr_comparison_plot()` figure per chromosome in a `for` loop. The new function is added to the existing `R/functions.R`, following the established pattern of the other plot functions in that file.

**Tech Stack:** R 4.5, Quarto, plotly, dplyr, purrr, readr, here, testthat 3, pixi

---

## File Map

| File | Action | Responsibility |
|---|---|---|
| `tests/test_make_chr_comparison_plot.R` | Create | Unit tests for the new plot function |
| `R/functions.R` | Modify | Add `make_chr_comparison_plot()` |
| `params_comparison.yml` | Create | Quarto params for both design directories |
| `chr_comparison_report.qmd` | Create | Report: setup chunk + plots loop |
| `pixi.toml` | Modify | Add `render-comparison` task |

---

## Task 1: Write the failing test for `make_chr_comparison_plot`

**Files:**
- Create: `tests/test_make_chr_comparison_plot.R`

- [ ] **Step 1: Create the test file**

```r
# tests/test_make_chr_comparison_plot.R
library(testthat)
library(tibble)
library(dplyr)
library(plotly)

source(file.path("..", "R", "functions.R"))

# ---------------------------------------------------------------------------
# Shared test fixture
# ---------------------------------------------------------------------------

make_test_bed <- function() {
  tibble(
    X1     = c("chr1", "chr1", "chr1", "chr1", "chr2", "chr2"),
    X2     = c(100L,   200L,   100L,   200L,   100L,   200L),
    X3     = c(200L,   300L,   200L,   300L,   200L,   300L),
    X4     = c("r1",   "r2",   "r1",   "r2",   "r3",   "r4"),
    X5     = c(50,     60,     70,     80,     90,     100),
    sample = c("S1",   "S1",   "S2",   "S2",   "S1",   "S2"),
    design = c("v1",   "v1",   "v2",   "v2",   "v1",   "v2")
  )
}

# ---------------------------------------------------------------------------
# 1. Happy path
# ---------------------------------------------------------------------------

test_that("make_chr_comparison_plot returns a plotly object", {
  p <- make_chr_comparison_plot("chr1", make_test_bed())
  expect_s3_class(p, "plotly")
})

# ---------------------------------------------------------------------------
# 2. Plot title contains the chromosome name
# ---------------------------------------------------------------------------

test_that("plot title contains the chromosome name", {
  result <- plotly_build(make_chr_comparison_plot("chr1", make_test_bed()))
  expect_match(result$x$layout$title, "chr1")
})

# ---------------------------------------------------------------------------
# 3. Warns when chromosome matches no rows
# ---------------------------------------------------------------------------

test_that("warns when chromosome matches no rows", {
  expect_warning(
    make_chr_comparison_plot("chr99", make_test_bed()),
    regexp = "no data found for chromosome"
  )
})

# ---------------------------------------------------------------------------
# 4. Invalid chromosome argument
# ---------------------------------------------------------------------------

test_that("errors when chromosome is not a single non-empty string", {
  bed <- make_test_bed()
  expect_error(make_chr_comparison_plot("",              bed), "single non-empty character string")
  expect_error(make_chr_comparison_plot(1L,              bed), "single non-empty character string")
  expect_error(make_chr_comparison_plot(c("chr1","chr2"),bed), "single non-empty character string")
  expect_error(make_chr_comparison_plot(NA_character_,   bed), "single non-empty character string")
})

# ---------------------------------------------------------------------------
# 5. Missing required columns in bed
# ---------------------------------------------------------------------------

test_that("errors when bed is missing required columns", {
  expect_error(
    make_chr_comparison_plot("chr1", make_test_bed() |> select(-X5)),
    regexp = "'bed' is missing required columns:.*X5"
  )
  expect_error(
    make_chr_comparison_plot("chr1", make_test_bed() |> select(-design)),
    regexp = "'bed' is missing required columns:.*design"
  )
  expect_error(
    make_chr_comparison_plot("chr1", make_test_bed() |> select(-X1)),
    regexp = "'bed' is missing required columns:.*X1"
  )
})
```

- [ ] **Step 2: Run the test to confirm it fails**

```bash
pixi run Rscript -e "testthat::test_file('tests/test_make_chr_comparison_plot.R')"
```

Expected: errors because `make_chr_comparison_plot` does not exist yet.

---

## Task 2: Implement `make_chr_comparison_plot` in `R/functions.R`

**Files:**
- Modify: `R/functions.R` (append after the last function, before end of file)

- [ ] **Step 1: Add the function to `R/functions.R`**

Add this block at the end of `R/functions.R` (after `join_bed_run_ids`):

```r
#' Plot Chromosome Coverage Comparison
#'
#' Generates an interactive Plotly boxplot comparing depth of coverage
#' between two experiment designs for a given chromosome.
#'
#' @param chromosome Character string. The chromosome name to filter for (e.g., "chr1").
#' @param bed Data frame. Combined BED data for both designs, with a 'design' column.
#'
#' @return A plotly object representing the boxplot.
#'
#' @export
#'
#' @examples
#' make_chr_comparison_plot("chr1", bed)
make_chr_comparison_plot <- function(chromosome, bed) {
  if (!is.character(chromosome) || length(chromosome) != 1 ||
      is.na(chromosome) || nchar(trimws(chromosome)) == 0) {
    stop("'chromosome' must be a single non-empty character string.")
  }

  required_cols <- c("X1", "X5", "design")
  missing_cols <- setdiff(required_cols, colnames(bed))
  if (length(missing_cols) > 0) {
    stop(
      "'bed' is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  chr_bed <- bed |> filter(.data$X1 == chromosome)

  if (nrow(chr_bed) == 0) {
    warning("make_chr_comparison_plot: no data found for chromosome '", chromosome, "'.")
  }

  plot_ly(
    data = chr_bed,
    x = ~design,
    y = ~X5,
    color = ~design,
    type = "box",
    boxpoints = "all",
    jitter = 0.3,
    pointpos = 0,
    marker = list(size = 4, opacity = 0.6),
    line = list(color = "rgba(0,0,0,0.5)")
  ) |>
    layout(
      title = paste0(chromosome, " depth of coverage"),
      yaxis = list(title = "depth"),
      xaxis = list(title = "design")
    )
}
```

- [ ] **Step 2: Run the tests to confirm they all pass**

```bash
pixi run Rscript -e "testthat::test_file('tests/test_make_chr_comparison_plot.R')"
```

Expected: all tests PASS, no failures.

- [ ] **Step 3: Run the full test suite to confirm no regressions**

```bash
pixi run test
```

Expected: all existing tests still PASS.

- [ ] **Step 4: Commit**

```bash
git add R/functions.R tests/test_make_chr_comparison_plot.R
git commit -m "feat: add make_chr_comparison_plot and its tests"
```

---

## Task 3: Create `params_comparison.yml`

**Files:**
- Create: `params_comparison.yml`

- [ ] **Step 1: Create the file**

```yaml
bed_dir_v1: "/path/to/design1/mosdepth_output"
bed_dir_v2: "/path/to/design2/mosdepth_output"
```

Replace the placeholder paths with the actual directories containing `.regions.bed.gz` files for each design before rendering.

- [ ] **Step 2: Commit**

```bash
git add params_comparison.yml
git commit -m "chore: add params_comparison.yml with placeholder paths"
```

---

## Task 4: Create `chr_comparison_report.qmd`

**Files:**
- Create: `chr_comparison_report.qmd`

- [ ] **Step 1: Create the report file**

````markdown
---
title: "Chromosome Coverage Comparison"
subtitle: "Design v1 vs v2"
author: "A.G."
format:
  html:
    self-contained: true
    theme: "cosmo"
    toc: true
    toc-location: left
editor: source
params:
  bed_dir_v1: "/path/to/v1"
  bed_dir_v2: "/path/to/v2"
---

```{r}
#| label: setup
#| echo: false
#| message: false
#| warning: false
#| cache: false
#| include: false

library(here)
library(dplyr)
library(purrr)
library(readr)
library(htmltools)

source(here("R", "functions.R"))

read_design <- function(dir, label) {
  list_bed_files(dir, ".regions.bed.gz") |>
    map_dfr(read_bed) |>
    mutate(design = label)
}

bed <- bind_rows(
  read_design(params$bed_dir_v1, "v1"),
  read_design(params$bed_dir_v2, "v2")
)

chromosomes <- sort(unique(bed$X1))
```

```{r}
#| label: plots
#| echo: false
#| results: "asis"

for (chr in chromosomes) {
  cat(paste0("\n## ", chr, "\n\n"))
  print(tagList(make_chr_comparison_plot(chr, bed)))
  cat("\n\n")
}
```
````

- [ ] **Step 2: Verify the report renders (requires real paths in params_comparison.yml)**

```bash
pixi run render-comparison
```

Expected: an HTML file named `coverage_comparison_<timestamp>.html` is created in the project root, with one section per chromosome, each containing a two-box Plotly figure.

- [ ] **Step 3: Commit**

```bash
git add chr_comparison_report.qmd
git commit -m "feat: add chromosome comparison Quarto report"
```

---

## Task 5: Add `render-comparison` task to `pixi.toml`

**Files:**
- Modify: `pixi.toml`

- [ ] **Step 1: Add the task**

In `pixi.toml`, under `[tasks]`, add:

```toml
render-comparison = "quarto render chr_comparison_report.qmd --output \"coverage_comparison_$(date +%Y-%m-%d_%H-%M-%S).html\" --execute-params params_comparison.yml"
```

The `[tasks]` block after the change should look like:

```toml
[tasks]
render-v2          = "quarto render main_report.qmd --output \"coverage_report_v2_$(date +%Y-%m-%d_%H-%M-%S).html\" --execute-params params_v2.yml"
render-v1          = "quarto render main_report.qmd --output \"coverage_report_v1_$(date +%Y-%m-%d_%H-%M-%S).html\" --execute-params params_v1.yml"
render-comparison  = "quarto render chr_comparison_report.qmd --output \"coverage_comparison_$(date +%Y-%m-%d_%H-%M-%S).html\" --execute-params params_comparison.yml"
test               = "Rscript -e 'testthat::test_dir(\"tests\")'"
style              = "Rscript -e 'styler::style_dir(\"R\")' -e 'styler::style_file(\"main_report.qmd\")' -e 'styler::style_dir(\"tests\")'"
lint               = "Rscript -e 'if(length(lintr::lint_dir(\"R\")) > 0) stop(\"Linting errors found!\")'"
```

- [ ] **Step 2: Commit**

```bash
git add pixi.toml
git commit -m "chore: add render-comparison pixi task"
```

---

## Self-Review

**Spec coverage:**
- [x] Two dynamic params (`bed_dir_v1`, `bed_dir_v2`) → Task 3
- [x] Directories of `.regions.bed.gz` files per design → Task 4 setup chunk
- [x] Combined Plotly boxplot per chromosome, x = design, y = coverage → Task 2 + Task 4
- [x] Shared y-axis (single figure with two traces) → satisfied by single `plot_ly()` call with `color = ~design`
- [x] `make_chr_comparison_plot` in `R/functions.R` → Task 2
- [x] Tests with same structure as existing test files → Task 1
- [x] `params_comparison.yml` → Task 3
- [x] `chr_comparison_report.qmd` → Task 4
- [x] `render-comparison` pixi task → Task 5

**No placeholders found.**

**Type consistency:** `make_chr_comparison_plot(chromosome, bed)` is used consistently in Task 1 (tests), Task 2 (implementation), and Task 4 (report).
