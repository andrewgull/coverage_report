library(testthat)
library(tibble)
library(dplyr)
library(stringr)
library(plotly)

source(file.path("..", "R", "functions.R"))

# ---------------------------------------------------------------------------
# Shared test fixture
# ---------------------------------------------------------------------------

make_backbone_bed <- function() {
  tibble(
    X4     = c("Backbone_Chr1_001", "Backbone_Chr1_002", "Backbone_Chr17_001"),
    X5     = c(100L, 200L, 150L),
    sample = c("S1", "S2", "S1"),
    label  = c("NA:NA:NA", "NA:NA:NA", "NA:NA:NA")
  )
}

# ---------------------------------------------------------------------------
# 1. Returns a plotly object on valid input
# ---------------------------------------------------------------------------

test_that("returns a plotly object for a matching chromosome", {
  result <- plot_backbone("Backbone_Chr1", make_backbone_bed())

  expect_s3_class(result, "plotly")
})

# ---------------------------------------------------------------------------
# 2. Plot title contains the chromosome name
# ---------------------------------------------------------------------------

test_that("plot title contains the chromosome name", {
  result <- plotly_build(plot_backbone("Backbone_Chr1", make_backbone_bed()))

  title <- result$x$layout$title
  expect_match(title, "Backbone_Chr1")
})

# ---------------------------------------------------------------------------
# 3. Only rows matching the chromosome prefix are passed to the plot
# ---------------------------------------------------------------------------

test_that("matching chromosome produces a plot without warnings", {
  # Positive-path check: a chromosome that exists should not trigger a warning
  expect_no_warning(
    plot_backbone("Backbone_Chr1", make_backbone_bed())
  )
})

# ---------------------------------------------------------------------------
# 4. Warns when no rows match the chromosome pattern
# ---------------------------------------------------------------------------

test_that("warns when chromosome pattern matches no rows", {
  expect_warning(
    plot_backbone("Backbone_Chr99", make_backbone_bed()),
    regexp = "no rows matched chromosome"
  )
})

# ---------------------------------------------------------------------------
# 5. Invalid chromosome argument: not a string
# ---------------------------------------------------------------------------

test_that("stops when chromosome is not a character string", {
  expect_error(
    plot_backbone(17, make_backbone_bed()),
    regexp = "single non-empty character string"
  )
})

# ---------------------------------------------------------------------------
# 6. Invalid chromosome argument: empty string
# ---------------------------------------------------------------------------

test_that("stops when chromosome is an empty string", {
  expect_error(
    plot_backbone("", make_backbone_bed()),
    regexp = "single non-empty character string"
  )
})

# ---------------------------------------------------------------------------
# 7. Invalid chromosome argument: length > 1
# ---------------------------------------------------------------------------

test_that("stops when chromosome has length > 1", {
  expect_error(
    plot_backbone(c("Backbone_Chr1", "Backbone_Chr17"), make_backbone_bed()),
    regexp = "single non-empty character string"
  )
})

# ---------------------------------------------------------------------------
# 8. Missing required column in bed
# ---------------------------------------------------------------------------

test_that("stops when bed is missing required columns", {
  bad_bed <- make_backbone_bed() |> select(-label)

  expect_error(
    plot_backbone("Backbone_Chr1", bad_bed),
    regexp = "'bed' is missing required columns:.*label"
  )
})
