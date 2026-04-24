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
    regexp = "make_chr_comparison_plot: no data found for chromosome"
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
  expect_error(make_chr_comparison_plot(NULL,            bed), "single non-empty character string")
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
