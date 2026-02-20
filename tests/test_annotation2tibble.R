library(testthat)
library(tibble)
library(stringr)

source(file.path("..", "R", "functions.R"))

# ---------------------------------------------------------------------------
# 1. Happy path: valid annotation string
# ---------------------------------------------------------------------------

test_that("annotation2tibble parses valid annotation strings correctly", {
    ann <- 'gene_name "TP53"; exon_number 1; transcript_name "TP53-201"; db_xref "RefSeq:NM_000546";'
    result <- annotation2tibble(ann)

    expect_s3_class(result, "tbl_df")
    expect_equal(result$gene_name, "TP53")
    expect_equal(result$exon_number, 1L)
    expect_equal(result$transcript_name, "TP53-201")
    expect_equal(result$refseq, "NM_000546")
})

# ---------------------------------------------------------------------------
# 2. Vector input
# ---------------------------------------------------------------------------

test_that("annotation2tibble handles vector of strings", {
    anns <- c(
        'gene_name "TP53"; exon_number 1; transcript_name "T1"; db_xref "RefSeq:R1";',
        'gene_name "KRAS"; exon_number 2; transcript_name "T2"; db_xref "RefSeq:R2";'
    )
    result <- annotation2tibble(anns)

    expect_equal(nrow(result), 2L)
    expect_equal(result$gene_name, c("TP53", "KRAS"))
    expect_equal(result$exon_number, c(1L, 2L))
})

# ---------------------------------------------------------------------------
# 3. Missing fields (regex doesn't match)
# ---------------------------------------------------------------------------

test_that("annotation2tibble returns NA for missing fields", {
    ann <- 'gene_name "TP53";' # missing others
    result <- annotation2tibble(ann)

    expect_equal(result$gene_name, "TP53")
    expect_true(is.na(result$exon_number))
    expect_true(is.na(result$transcript_name))
    expect_true(is.na(result$refseq))
})

# ---------------------------------------------------------------------------
# 4. Invalid inputs
# ---------------------------------------------------------------------------

test_that("annotation2tibble throws error for invalid inputs", {
    expect_error(annotation2tibble(NULL), "'column' must not be NULL.")
    expect_error(annotation2tibble(123), "'column' must be a character vector.")
})

# ---------------------------------------------------------------------------
# 5. Empty vector
# ---------------------------------------------------------------------------

test_that("annotation2tibble returns empty tibble for empty vector", {
    result <- annotation2tibble(character())
    expect_equal(nrow(result), 0L)
    expect_named(result, c("gene_name", "exon_number", "transcript_name", "refseq"))
})
