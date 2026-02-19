library(tidyr)
library(dplyr)
library(readr)
library(ggplot2)
library(purrr)
library(DT)
library(plotly)
library(stringr)
library(tibble)


annotation2tibble <- function(column) {
  # transform the targets
  # turn column with annotations (X13) into tibble
  # column: vector of character strings
  
  # extract gene name
  gene_name <- str_extract(column, 'gene_name "([^"]+)";') |>
    str_remove('gene_name "') |>
    str_remove('";')
  # extract exon number
  exon_number <- as.integer(
    str_extract(column, 'exon_number ([^"]+);') |>
      str_remove('exon_number ') |>
      str_remove(';')
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
  
  output <- tibble(
    "gene_name" = gene_name,
    "exon_number" = exon_number,
    "transcript_name" = transcript_name,
    "refseq" = refseq
  )
  
  return(output)
}


read_bed <- function(filename) {
  # read a BED file
  # filename: name of a BED file
  
  # extract sample name first
  # replace underscores
  sample_name <- basename(filename) |>
    str_remove("\\.regions\\.bed\\.gz$")

  # read bed file and add the sample name
  bed <- read_tsv(filename, col_names = F, show_col_types = F) |>
    mutate("sample" = sample_name)
  if (nrow(bed) == 0) {
    print(paste0("empty file ", filename))
  }
  
  return(bed)
}


make_gene_plot <- function(gene, bed = bed_annotated) {
  # gene: gene name (character string)
  # bed: bed file object (data frame/tibble)
  
  bed_gene <- bed |> dplyr::filter(gene_name == gene)
  
  fig <- plot_ly(
    data = bed_gene,
    y = ~ X5,
    x = ~ label,
    type = "box",
    boxpoints = "all",
    jitter = 0.3,
    pointpos = 0,
    marker = list(size = 4, opacity = 0.6),
    line = list(color = 'rgba(0,0,0,0.5)'),
    fillcolor = 'transparent',
    text = ~ sample
  ) |>
    layout(
      title = paste0(gene, " median depth of coverage"),
      yaxis = list(title = "depth"),
      xaxis = list(title = "target region")
    )
  
  return(fig)
}


make_sample_plot <- function(bed,
                             group_by_gene = FALSE,
                             gene = "") {
  # bed: bed file object (data frame/tibble)
  # group_by_gene: group by gene or not (logical)
  # gene: gene to make plot for (character string, only if group_by_gene=TRUE)
  if (!group_by_gene) {
    plot_ly(
      data = bed,
      y = ~ X5,
      x = ~ sample,
      type = "box",
      boxpoints = "all",
      jitter = 0.3,
      pointpos = 0,
      marker = list(
        size = 4,
        opacity = 0.6,
        color = 'rgba(102,0,0,0.5)'
      ),
      line = list(color = 'rgba(102,0,0,0.5)'),
      fillcolor = 'transparent',
      text = ~ label
    ) |>
      layout(
        title = "Median depth of coverage",
        yaxis = list(title = "depth"),
        xaxis = list(title = "sample")
      )
  } else {
    plot_ly(
      data = dplyr::filter(bed, gene_name == gene),
      y = ~ X5,
      x = ~ sample,
      type = "box",
      boxpoints = "all",
      jitter = 0.3,
      pointpos = 0,
      marker = list(
        size = 4,
        opacity = 0.6,
        color = 'rgba(102,0,0,0.5)'
      ),
      line = list(color = 'rgba(102,0,0,0.5)'),
      fillcolor = 'transparent',
      text = ~ label
    ) |>
      layout(
        title = paste0(gene, " median depth of coverage"),
        yaxis = list(title = "depth"),
        xaxis = list(title = "sample")
      )
  }
  
}


plot_backbone <- function(chromosome, bed = bed_annotated_backbone) {
  # chromosome: chromosome name (character sting: chr1 or ch17)
  # bed: bed object with backbone targets (data frame)
  chr_bed <- bed |>
    filter(grepl(str_glue("^{chromosome}_"), X4))
  
  fig <- plot_ly(
    data = chr_bed,
    y = ~ X5,
    x = ~ sample,
    type = "box",
    boxpoints = "all",
    jitter = 0.3,
    pointpos = 0,
    marker = list(
      size = 4,
      opacity = 0.6,
      color = 'rgba(102,0,0,0.5)'
    ),
    line = list(color = 'rgba(102,0,0,0.5)'),
    fillcolor = 'transparent',
    text = ~ label
  ) |>
    layout(
      title = paste0(chromosome, " median depth of coverage"),
      yaxis = list(title = "depth"),
      xaxis = list(title = "sample")
    )
  
  return(fig)
}


make_cov_table <- function(bed, by_gene = FALSE) {
  # bed: data frame object representing BED file
  if (by_gene) {
    bed |>
      group_by(gene_name, label) |>
      summarise(
        median_coverage = round(median(X5), 0),
        min_coverage = round(min(X5), 0),
        max_coverage = round(max(X5), 0)
      )
  } else {
    bed |>
      group_by(sample) |>
      summarise(
        median_coverage = round(median(X5), 0),
        min_coverage = round(min(X5), 0),
        max_coverage = round(max(X5), 0)
      )
  }
  
}
