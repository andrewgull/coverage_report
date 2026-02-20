library(testthat)
library(tibble)
library(readr)
library(dplyr)
library(stringr)

source(file.path("..", "R", "functions.R"))

# ---------------------------------------------------------------------------
# 1. Happy path: valid BED file
# ---------------------------------------------------------------------------

test_that("read_bed returns a data frame with correct sample name", {
  tmp_file <- tempfile(fileext = ".regions.bed.gz")
  test_data <- tibble(X1 = "chr1", X2 = 100, X3 = 200, X4 = "Target", X5 = 50)
  write_tsv(test_data, tmp_file, col_names = FALSE)

  result <- read_bed(tmp_file)

  expected_sample <- basename(tmp_file) |> str_remove("\\.regions\\.bed\\.gz$")
  expect_s3_class(result, "tbl_df")
  expect_equal(result$sample[1], expected_sample)
  expect_equal(nrow(result), 1L)

  unlink(tmp_file)
})

# ---------------------------------------------------------------------------
# 2. Invalid filename parameter
# ---------------------------------------------------------------------------

test_that("read_bed throws error for invalid filename", {
  expect_error(read_bed(123), "'filename' must be a single non-empty character string.")
  expect_error(read_bed(""), "'filename' must be a single non-empty character string.")
  expect_error(read_bed(c("file1", "file2")), "'filename' must be a single non-empty character string.")
})

# ---------------------------------------------------------------------------
# 3. File not found
# ---------------------------------------------------------------------------

test_that("read_bed throws error when file does not exist", {
  expect_error(read_bed("non_existent_file.regions.bed.gz"), "File '.*' does not exist.")
})

# ---------------------------------------------------------------------------
# 4. Empty file
# ---------------------------------------------------------------------------

test_that("read_bed warns when file is empty", {
  tmp_file <- tempfile(fileext = ".regions.bed.gz")
  file.create(tmp_file)

  expect_warning(read_bed(tmp_file), "read_bed: file '.*' is empty.")

  unlink(tmp_file)
})
