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
    gene_name = c("G1", "G1", "G2", "G2"),
    X5        = c(10, 20, 15, 25),
    label     = c("L1", "L2", "L1", "L2"),
    sample    = c("S1", "S1", "S2", "S2")
  )
}

# ---------------------------------------------------------------------------
# 1. Happy path
# ---------------------------------------------------------------------------

test_that("make_gene_plot returns a plotly object", {
  bed <- make_test_bed()
  p <- make_gene_plot("G1", bed)
  expect_s3_class(p, "plotly")
})

# ---------------------------------------------------------------------------
# 2. Invalid gene parameter
# ---------------------------------------------------------------------------

test_that("make_gene_plot throws error for invalid gene", {
  bed <- make_test_bed()
  expect_error(make_gene_plot(123, bed), "'gene' must be a single non-empty character string.")
  expect_error(make_gene_plot("", bed), "'gene' must be a single non-empty character string.")
  expect_error(make_gene_plot(c("G1", "G2"), bed), "'gene' must be a single non-empty character string.")
})

# ---------------------------------------------------------------------------
# 3. Missing columns in bed
# ---------------------------------------------------------------------------

test_that("make_gene_plot throws error for missing columns", {
  bed <- make_test_bed() |> select(-X5)
  expect_error(make_gene_plot("G1", bed), "'bed' is missing required columns: X5")
})

# ---------------------------------------------------------------------------
# 4. Gene not found / Warning
# ---------------------------------------------------------------------------

test_that("make_gene_plot warns when gene is not found", {
  bed <- make_test_bed()
  expect_warning(make_gene_plot("G3", bed), "make_gene_plot: no data found for gene 'G3'.")
})
