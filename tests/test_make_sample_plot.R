library(testthat)
library(tibble)
library(dplyr)
library(plotly)

source(file.path("..", "R", "functions.R"))

# ---------------------------------------------------------------------------
# Shared test fixtures
# ---------------------------------------------------------------------------

make_test_bed <- function() {
  tibble(
    sample    = c("S1", "S1", "S2", "S2"),
    X5        = c(10, 20, 15, 25),
    label     = c("L1", "L2", "L1", "L2"),
    gene_name = c("G1", "G1", "G2", "G2")
  )
}

# ---------------------------------------------------------------------------
# 1. Happy path: by sample
# ---------------------------------------------------------------------------

test_that("make_sample_plot returns a plotly object (by sample)", {
  bed <- make_test_bed()
  p <- make_sample_plot(bed, group_by_gene = FALSE)
  expect_s3_class(p, "plotly")
})

# ---------------------------------------------------------------------------
# 2. Happy path: by gene
# ---------------------------------------------------------------------------

test_that("make_sample_plot returns a plotly object (by gene)", {
  bed <- make_test_bed()
  p <- make_sample_plot(bed, group_by_gene = TRUE, gene = "G1")
  expect_s3_class(p, "plotly")
})

# ---------------------------------------------------------------------------
# 3. Invalid parameter types
# ---------------------------------------------------------------------------

test_that("make_sample_plot throws error for invalid group_by_gene", {
  bed <- make_test_bed()
  expect_error(make_sample_plot(bed, group_by_gene = "TRUE"), "'group_by_gene' must be logical.")
})

test_that("make_sample_plot throws error for invalid gene when group_by_gene is TRUE", {
  bed <- make_test_bed()
  expect_error(make_sample_plot(bed, group_by_gene = TRUE, gene = ""), "When 'group_by_gene' is TRUE, 'gene' must be a non-empty character string.")
  expect_error(make_sample_plot(bed, group_by_gene = TRUE, gene = 123), "When 'group_by_gene' is TRUE, 'gene' must be a non-empty character string.")
})

# ---------------------------------------------------------------------------
# 4. Missing columns
# ---------------------------------------------------------------------------

test_that("make_sample_plot throws error for missing columns (by sample)", {
  bed <- make_test_bed() |> select(-X5)
  expect_error(make_sample_plot(bed, group_by_gene = FALSE), "'bed' is missing required columns: X5")
})

test_that("make_sample_plot throws error for missing columns (by gene)", {
  bed <- make_test_bed() |> select(-gene_name)
  expect_error(make_sample_plot(bed, group_by_gene = TRUE, gene = "G1"), "'bed' is missing required columns: gene_name")
})

# ---------------------------------------------------------------------------
# 5. Empty data / Warnings
# ---------------------------------------------------------------------------

test_that("make_sample_plot warns when bed is empty (by sample)", {
  bed <- make_test_bed() |> filter(sample == "NONE")
  expect_warning(make_sample_plot(bed, group_by_gene = FALSE), "make_sample_plot: input 'bed' is empty.")
})

test_that("make_sample_plot warns when gene is not found (by gene)", {
  bed <- make_test_bed()
  expect_warning(make_sample_plot(bed, group_by_gene = TRUE, gene = "NONEXISTENT"), "make_sample_plot: no data found for gene 'NONEXISTENT'.")
})
