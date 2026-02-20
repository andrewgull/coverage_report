library(tidyr)
library(dplyr)
library(readr)
library(ggplot2)
library(purrr)
library(DT)
library(plotly)
library(stringr)
library(tibble)
library(rlang)

# Declare global variables to avoid 'no visible binding' warnings
utils::globalVariables(c(
  "gene_name", "sample",
  "label", "X4", "X5", "exon_number", "refseq",
  "run_ID", "transcript_name"
))


#' Parse annotation column
#'
#' @param column Character string. The annotation column to parse.
#' @return A data frame containing the parsed annotation.
#' @export
#' @examples
#' annotation2tibble("gene_name \"TP53\"; exon_number 1; transcript_name \"TP53\"; db_xref \"RefSeq:12345\";")
annotation2tibble <- function(column) {
  if (is.null(column)) {
    stop("'column' must not be NULL.")
  }
  if (!is.character(column)) {
    stop("'column' must be a character vector.")
  }

  if (length(column) == 0) {
    return(tibble(
      "gene_name" = character(),
      "exon_number" = integer(),
      "transcript_name" = character(),
      "refseq" = character()
    ))
  }

  # extract gene name
  gene_name <- str_extract(column, 'gene_name "([^"]+)";') |>
    str_remove('gene_name "') |>
    str_remove('";')
  # extract exon number
  exon_number <- as.integer(
    str_extract(column, 'exon_number ([^"]+);') |>
      str_remove("exon_number ") |>
      str_remove(";")
  )
  # extract transcript name
  transcript_name <- str_extract(column, 'transcript_name "([^"]+)";') |>
    str_remove('transcript_name "') |>
    str_remove('";')
  # extract refseq
  refseq <- str_extract(column, 'db_xref "([^"]+)";') |>
    str_remove('db_xref "') |>
    str_remove('";') |>
    str_remove("RefSeq:")

  tibble(
    "gene_name" = gene_name,
    "exon_number" = exon_number,
    "transcript_name" = transcript_name,
    "refseq" = refseq
  )
}

#' Parse BED file
#'
#' @param filename Character string. The name of the BED file to parse.
#' @return A data frame containing the parsed BED file.
#' @export
#' @examples
#' read_bed("sample1.regions.bed.gz")
read_bed <- function(filename) {
  if (!is.character(filename) || length(filename) != 1 || is.na(filename) || nchar(trimws(filename)) == 0) {
    stop("'filename' must be a single non-empty character string.")
  }

  if (!file.exists(filename)) {
    stop("File '", filename, "' does not exist.")
  }

  # extract sample name first
  # replace underscores
  sample_name <- basename(filename) |>
    str_remove("\\.regions\\.bed\\.gz$")

  # read bed file and add the sample name
  bed <- tryCatch(
    {
      read_tsv(filename, col_names = FALSE, show_col_types = FALSE)
    },
    error = function(e) {
      stop("Error reading file '", filename, "': ", e$message)
    }
  )

  if (nrow(bed) == 0) {
    warning("read_bed: file '", filename, "' is empty.")
  }

  bed |>
    mutate("sample" = sample_name)
}

#' Plot Gene Coverage
#'
#' Generates an interactive Plotly boxplot showing the median depth of coverage
#' for a specific gene based on annotated targets.
#'
#' @param gene Character string. The gene name to filter for (e.g., "TP53").
#' @param bed Data frame. A BED-formatted object containing annotated targets.
#'   Defaults to `bed_annotated`.
#'
#' @return A plotly object representing the boxplot.
#'
#' @export
#'
#' @examples
#' plot_gene("TP53")
make_gene_plot <- function(gene, bed = bed_annotated) {
  if (!is.character(gene) || length(gene) != 1 || nchar(trimws(gene)) == 0) {
    stop("'gene' must be a single non-empty character string.")
  }

  required_cols <- c("gene_name", "X5", "label", "sample")
  missing_cols <- setdiff(required_cols, colnames(bed))
  if (length(missing_cols) > 0) {
    stop(
      "'bed' is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  bed_gene <- bed |> dplyr::filter(.data$gene_name == gene)

  if (nrow(bed_gene) == 0) {
    warning("make_gene_plot: no data found for gene '", gene, "'.")
  }

  fig <- plot_ly(
    data = bed_gene,
    y = ~X5,
    x = ~label,
    type = "box",
    boxpoints = "all",
    jitter = 0.3,
    pointpos = 0,
    marker = list(size = 4, opacity = 0.6),
    line = list(color = "rgba(0,0,0,0.5)"),
    fillcolor = "transparent",
    text = ~sample
  ) |>
    layout(
      title = paste0(gene, " median depth of coverage"),
      yaxis = list(title = "depth"),
      xaxis = list(title = "target region")
    )

  fig
}

#' Plot Sample Coverage
#'
#' Generates an interactive Plotly boxplot showing the median depth of coverage
#' for a specific sample based on annotated backbone targets.
#'
#' @param bed Data frame. A BED-formatted object containing backbone targets.
#'   Defaults to `bed_annotated_backbone`.
#' @param group_by_gene Logical. Whether to group by gene or not.
#' @param gene Character string. The gene to make plot for (only if group_by_gene=TRUE).
#'
#' @return A plotly object representing the boxplot.
#'
#' @export
#'
#' @examples
#' plot_sample_coverage()
make_sample_plot <- function(bed,
                             group_by_gene = FALSE,
                             gene = "") {
  # bed: bed file object (data frame/tibble)
  # group_by_gene: group by gene or not (logical)
  # gene: gene to make plot for (character string, only if group_by_gene=TRUE)

  if (!is.logical(group_by_gene)) {
    stop("'group_by_gene' must be logical.")
  }

  required_cols <- c("X5", "sample", "label")
  if (group_by_gene) {
    required_cols <- c(required_cols, "gene_name")
    if (!is.character(gene) || length(gene) != 1 || is.na(gene) || nchar(trimws(gene)) == 0) {
      stop("When 'group_by_gene' is TRUE, 'gene' must be a non-empty character string.")
    }
  }

  missing_cols <- setdiff(required_cols, colnames(bed))
  if (length(missing_cols) > 0) {
    stop(
      "'bed' is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  plot_data <- if (group_by_gene) {
    dplyr::filter(bed, .data$gene_name == gene)
  } else {
    bed
  }

  if (nrow(plot_data) == 0) {
    msg <- if (group_by_gene) {
      paste0("make_sample_plot: no data found for gene '", gene, "'.")
    } else {
      "make_sample_plot: input 'bed' is empty."
    }
    warning(msg)
  }

  if (!group_by_gene) {
    plot_ly(
      data = plot_data,
      y = ~X5,
      x = ~sample,
      type = "box",
      boxpoints = "all",
      jitter = 0.3,
      pointpos = 0,
      marker = list(
        size = 4,
        opacity = 0.6,
        color = "rgba(102,0,0,0.5)"
      ),
      line = list(color = "rgba(102,0,0,0.5)"),
      fillcolor = "transparent",
      text = ~label
    ) |>
      layout(
        title = "Median depth of coverage",
        yaxis = list(title = "depth"),
        xaxis = list(title = "sample")
      )
  } else {
    plot_ly(
      data = plot_data,
      y = ~X5,
      x = ~sample,
      type = "box",
      boxpoints = "all",
      jitter = 0.3,
      pointpos = 0,
      marker = list(
        size = 4,
        opacity = 0.6,
        color = "rgba(102,0,0,0.5)"
      ),
      line = list(color = "rgba(102,0,0,0.5)"),
      fillcolor = "transparent",
      text = ~label
    ) |>
      layout(
        title = paste0(gene, " median depth of coverage"),
        yaxis = list(title = "depth"),
        xaxis = list(title = "sample")
      )
  }
}

#' Plot Backbone Coverage
#'
#' Generates an interactive Plotly boxplot showing the median depth of coverage
#' for a specific chromosome based on annotated backbone targets.
#'
#' @param chromosome Character string. The chromosome name to filter for
#'   (e.g., "chr1" or "chr17").
#' @param bed Data frame. A BED-formatted object containing backbone targets.
#'   Defaults to `bed_annotated_backbone`.
#'
#' @return A plotly object representing the boxplot.
#'
#' @export
#'
#' @examples
#' plot_backbone("chr1")
plot_backbone <- function(chromosome, bed = bed_annotated_backbone) {
  if (!is.character(chromosome) || length(chromosome) != 1 ||
    nchar(trimws(chromosome)) == 0) { # nolint
    stop("'chromosome' must be a single non-empty character string.")
  }

  required_cols <- c("X4", "X5", "sample", "label")
  missing_cols <- setdiff(required_cols, colnames(bed))
  if (length(missing_cols) > 0) {
    stop(
      "'bed' is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  chr_bed <- bed |>
    filter(grepl(str_glue("^{chromosome}_"), .data$X4))

  if (nrow(chr_bed) == 0) {
    warning(
      "plot_backbone: no rows matched chromosome '", chromosome, "'. ",
      "Check that X4 values start with '", chromosome, "_'."
    )
  }

  fig <- plot_ly(
    data = chr_bed,
    y = ~X5,
    x = ~sample,
    type = "box",
    boxpoints = "all",
    jitter = 0.3,
    pointpos = 0,
    marker = list(
      size = 4,
      opacity = 0.6,
      color = "rgba(102,0,0,0.5)"
    ),
    line = list(color = "rgba(102,0,0,0.5)"),
    fillcolor = "transparent",
    text = ~label
  ) |>
    layout(
      title = paste0(chromosome, " median depth of coverage"),
      yaxis = list(title = "depth"),
      xaxis = list(title = "sample")
    )

  fig
}

#' Make Coverage Table
#'
#' Generates a table of coverage statistics from a BED-formatted object.
#'
#' @param bed Data frame. A BED-formatted object containing backbone targets.
#'   Defaults to `bed_annotated_backbone`.
#' @param by_gene Logical. Whether to group by gene or not.
#'
#' @return A data frame/tibble containing the coverage table.
#'
#' @export
#'
#' @examples
#' make_cov_table(bed_annotated)
make_cov_table <- function(bed, by_gene = FALSE) {
  required_cols <- if (by_gene) {
    c("X5", "gene_name", "label")
  } else {
    c("X5", "sample")
  }

  missing_cols <- setdiff(required_cols, colnames(bed))
  if (length(missing_cols) > 0) {
    stop(
      "'bed' is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  if (nrow(bed) == 0) {
    warning("make_cov_table: 'bed' has 0 rows — returning empty table.")
  }

  if (by_gene) {
    bed |>
      group_by(.data$gene_name, .data$label) |>
      summarise(
        median_coverage = round(median(.data$X5), 0),
        min_coverage = round(min(.data$X5), 0),
        max_coverage = round(max(.data$X5), 0),
        .groups = "drop"
      )
  } else {
    bed |>
      group_by(.data$sample) |>
      summarise(
        median_coverage = round(median(.data$X5), 0),
        min_coverage = round(min(.data$X5), 0),
        max_coverage = round(max(.data$X5), 0),
        .groups = "drop"
      )
  }
}

#' Join Annotated BED and RUN IDs Table
#'
#' @param bed Data frame. A BED-formatted object containing backbone targets.
#' @param run_ids Data frame. A data frame containing run IDs.
#'
#' @return A data frame/tibble containing the joined by sample data.
#' @export
#' @examples
#' join_bed_run_ids(bed_annotated, run_ids)
join_bed_run_ids <- function(bed, run_ids) {
  required_bed_cols <- c("sample", "X5", "transcript_name", "refseq")
  required_run_cols <- c("sample", "run_ID")

  missing_bed <- setdiff(required_bed_cols, colnames(bed))
  if (length(missing_bed) > 0) {
    stop(
      "'bed' is missing required columns: ",
      paste(missing_bed, collapse = ", ")
    )
  }

  missing_run <- setdiff(required_run_cols, colnames(run_ids))
  if (length(missing_run) > 0) {
    stop(
      "'run_ids' is missing required columns: ",
      paste(missing_run, collapse = ", ")
    )
  }

  if (length(intersect(bed$sample, run_ids$sample)) == 0) {
    stop(
      "No 'sample' values match between 'bed' and 'run_ids'. ",
      "Check for suffix mismatches (e.g. '_T' or '_N').\n",
      "  bed samples:     ",
      paste(head(unique(bed$sample), 5), collapse = ", "), "\n",
      "  run_ids samples: ",
      paste(head(unique(run_ids$sample), 5), collapse = ", ")
    )
  }

  result <- bed |>
    left_join(run_ids, by = "sample") |>
    select(X5, run_ID, sample, transcript_name, refseq)

  if (nrow(result) == 0) {
    warning(
      "join_bed_run_ids: result has 0 rows — ",
      "no 'sample' values matched between 'bed' and 'run_ids'."
    )
  }

  result
}
