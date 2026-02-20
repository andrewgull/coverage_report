library(testthat)
library(fs)

source(file.path("..", "R", "functions.R"))

# ---------------------------------------------------------------------------
# 1. Happy path: valid path and suffix with matching files
# ---------------------------------------------------------------------------

test_that("list_bed_files returns correct file paths when files exist", {
  tmp_dir <- tempfile("test_beds")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE))

  file1 <- file.path(tmp_dir, "sample1.regions.bed.gz")
  file2 <- file.path(tmp_dir, "sample2.regions.bed.gz")
  file.create(file1)
  file.create(file2)
  # A file that shouldn't match
  file.create(file.path(tmp_dir, "other.txt"))

  result <- list_bed_files(tmp_dir, ".regions.bed.gz")

  expect_type(result, "character")
  expect_length(result, 2)
  expect_true(all(basename(result) %in% c("sample1.regions.bed.gz", "sample2.regions.bed.gz")))
  expect_true(all(file.exists(result)))
})

# ---------------------------------------------------------------------------
# 2. Invalid 'path' parameter
# ---------------------------------------------------------------------------

test_that("list_bed_files throws error for invalid 'path'", {
  expect_error(list_bed_files(123, ".gz"), "'path' must be a single non-empty character string.")
  expect_error(list_bed_files(NA_character_, ".gz"), "'path' must be a single non-empty character string.")
  expect_error(list_bed_files("", ".gz"), "'path' must be a single non-empty character string.")
  expect_error(list_bed_files("  ", ".gz"), "'path' must be a single non-empty character string.")
  expect_error(list_bed_files(c("p1", "p2"), ".gz"), "'path' must be a single non-empty character string.")
})

# ---------------------------------------------------------------------------
# 3. Invalid 'suffix' parameter
# ---------------------------------------------------------------------------

test_that("list_bed_files throws error for invalid 'suffix'", {
  expect_error(list_bed_files(".", 123), "'suffix' must be a single non-empty character string.")
  expect_error(list_bed_files(".", NA_character_), "'suffix' must be a single non-empty character string.")
  expect_error(list_bed_files(".", ""), "'suffix' must be a single non-empty character string.")
  expect_error(list_bed_files(".", "  "), "'suffix' must be a single non-empty character string.")
  expect_error(list_bed_files(".", c("s1", "s2")), "'suffix' must be a single non-empty character string.")
})

# ---------------------------------------------------------------------------
# 4. No files found
# ---------------------------------------------------------------------------

test_that("list_bed_files throws error when no files match the suffix", {
  tmp_dir <- tempfile("test_no_match")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE))

  expect_error(list_bed_files(tmp_dir, ".regions.bed.gz"), "No files found with suffix")
})
