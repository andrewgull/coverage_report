library(testthat)
library(tibble)
library(dplyr)

source(file.path("..", "R", "functions.R"))

# ---------------------------------------------------------------------------
# Shared test fixtures
# ---------------------------------------------------------------------------

make_bed_sample <- function() {
  tibble(
    sample = c("S1", "S1", "S2", "S2"),
    X5     = c(100L, 200L, 50L, 150L)
  )
}

make_bed_gene <- function() {
  tibble(
    gene_name = c("TP53", "TP53", "KRAS", "KRAS"),
    label     = c("TP53:1:NM_1", "TP53:1:NM_1", "KRAS:2:NM_2", "KRAS:2:NM_2"),
    sample    = c("S1", "S2", "S1", "S2"),
    X5        = c(100L, 200L, 50L, 150L)
  )
}

# ---------------------------------------------------------------------------
# 1. by_gene = FALSE (default): groups by sample
# ---------------------------------------------------------------------------

test_that("default (by_gene=FALSE) groups by sample and returns correct stats", {
  result <- make_cov_table(make_bed_sample())

  expect_s3_class(result, "tbl_df")
  expect_named(
    result,
    c("sample", "median_coverage", "min_coverage", "max_coverage")
  )
  expect_equal(nrow(result), 2L)

  s1 <- result[result$sample == "S1", ]
  expect_equal(s1$median_coverage, 150L) # median(100, 200) = 150
  expect_equal(s1$min_coverage, 100L)
  expect_equal(s1$max_coverage, 200L)
})

# ---------------------------------------------------------------------------
# 2. by_gene = TRUE: groups by gene_name + label
# ---------------------------------------------------------------------------

test_that("by_gene=TRUE groups by gene_name and label and returns correct stats", {
  result <- make_cov_table(make_bed_gene(), by_gene = TRUE)

  expect_named(
    result,
    c("gene_name", "label", "median_coverage", "min_coverage", "max_coverage")
  )
  expect_equal(nrow(result), 2L) # two unique gene+label combos

  tp53 <- result[result$gene_name == "TP53", ]
  expect_equal(tp53$median_coverage, 150L) # median(100, 200)
})

# ---------------------------------------------------------------------------
# 3. Missing X5 column is always required
# ---------------------------------------------------------------------------

test_that("stops when X5 column is missing", {
  bad_bed <- tibble(sample = "S1")

  expect_error(
    make_cov_table(bad_bed),
    regexp = "'bed' is missing required columns:.*X5"
  )
})

# ---------------------------------------------------------------------------
# 4. Missing 'sample' when by_gene = FALSE
# ---------------------------------------------------------------------------

test_that("stops when 'sample' column missing and by_gene=FALSE", {
  bad_bed <- tibble(X5 = 100L)

  expect_error(
    make_cov_table(bad_bed, by_gene = FALSE),
    regexp = "'bed' is missing required columns:.*sample"
  )
})

# ---------------------------------------------------------------------------
# 5. Missing gene_name / label when by_gene = TRUE
# ---------------------------------------------------------------------------

test_that("stops when gene_name or label missing and by_gene=TRUE", {
  bad_bed <- tibble(X5 = 100L, gene_name = "TP53") # missing label

  expect_error(
    make_cov_table(bad_bed, by_gene = TRUE),
    regexp = "'bed' is missing required columns:.*label"
  )
})

# ---------------------------------------------------------------------------
# 6. Empty input produces a warning and an empty-row result
# ---------------------------------------------------------------------------

test_that("warns and returns 0-row tibble for empty input", {
  empty_bed <- make_bed_sample()[0, ]

  expect_warning(
    result <- make_cov_table(empty_bed),
    regexp = "0 rows"
  )
  expect_equal(nrow(result), 0L)
})

# ---------------------------------------------------------------------------
# 7. Coverage values are rounded to 0 decimal places
# ---------------------------------------------------------------------------

test_that("coverage values are integers (rounded to 0 decimals)", {
  bed <- tibble(sample = "S1", X5 = c(101L, 102L, 103L))
  result <- make_cov_table(bed)

  expect_equal(result$median_coverage, 102L)
  expect_true(result$median_coverage == as.integer(result$median_coverage))
})
