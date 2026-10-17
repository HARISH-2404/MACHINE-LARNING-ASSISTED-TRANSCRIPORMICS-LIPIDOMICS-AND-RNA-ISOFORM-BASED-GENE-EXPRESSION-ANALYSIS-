#!/usr/bin/env Rscript
# Isoform switching workflow using IsoformSwitchAnalyzeR (advanced; requires transcript-level estimates).
# This script is a template because input construction depends on the quantification and annotation.
# Usage: Rscript R/isoform_switch_analysis.R isoform_input.csv output_dir
# Required CSV columns: isoform_id,gene_id,sample_id,condition,abundance
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) stop("Usage: Rscript R/isoform_switch_analysis.R isoform_input.csv output_dir")
suppressPackageStartupMessages({library(readr); library(dplyr)})
if (!requireNamespace("IsoformSwitchAnalyzeR", quietly=TRUE)) {
  stop("Install IsoformSwitchAnalyzeR from Bioconductor before running this optional module.")
}
library(IsoformSwitchAnalyzeR)
d <- read_csv(args[1], show_col_types=FALSE); out <- args[2]; dir.create(out, recursive=TRUE, showWarnings=FALSE)
need <- c("isoform_id","gene_id","sample_id","condition","abundance")
if (!all(need %in% names(d))) stop("Input must contain: ", paste(need,collapse=", "))
if (anyNA(d[,need])) stop("Required columns contain missing values")
# Convert long-form abundances to matrices and verify one value per isoform/sample.
if (anyDuplicated(d[,c("isoform_id","sample_id")])) stop("Duplicate isoform_id/sample_id rows found")
counts <- tidyr::pivot_wider(d, id_cols=isoform_id, names_from=sample_id, values_from=abundance)
counts <- as.data.frame(counts); rownames(counts) <- counts$isoform_id; counts$isoform_id <- NULL
iso2gene <- d %>% distinct(isoform_id,gene_id)
write_csv(iso2gene, file.path(out,"isoform_to_gene_mapping.csv"))
write_csv(counts %>% tibble::rownames_to_column("isoform_id"), file.path(out,"isoform_abundance_matrix.csv"))
message("Prepared transcript abundance matrix and mapping. For a formal switch test, build a SwitchAnalyzeRlist with the matching transcript annotation and quantification method as documented by IsoformSwitchAnalyzeR.")
message("Do not interpret differential transcript abundance alone as proof of isoform switching.")