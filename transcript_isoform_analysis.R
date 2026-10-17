#!/usr/bin/env Rscript
# Simplified transcript-feature differential abundance analysis using integer transcript counts.
# This is NOT a complete isoform-switching or differential-splicing analysis.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript R/transcript_isoform_analysis.R transcript_counts.csv metadata.csv output_dir")
}
count_file <- args[1]; meta_file <- args[2]; out_dir <- args[3]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages({
  library(DESeq2); library(readr); library(dplyr); library(tibble)
})
df <- read_csv(count_file, show_col_types = FALSE)
meta <- read_csv(meta_file, show_col_types = FALSE)
if (!"feature_id" %in% names(df)) stop("Transcript matrix needs feature_id as first column.")
if (!all(c("sample_id", "condition") %in% names(meta))) stop("Metadata needs sample_id and condition.")
if (anyDuplicated(df$feature_id) || anyDuplicated(meta$sample_id)) stop("IDs must be unique.")
if (!setequal(setdiff(names(df), "feature_id"), meta$sample_id)) stop("Sample IDs do not match.")
df <- df %>% select(feature_id, all_of(meta$sample_id))
mat <- as.matrix(df[, -1]); rownames(mat) <- df$feature_id
storage.mode(mat) <- "numeric"
if (anyNA(mat) || any(mat < 0) || any(abs(mat - round(mat)) > 1e-8)) stop("Counts must be non-negative integers.")
storage.mode(mat) <- "integer"
meta <- as.data.frame(meta); rownames(meta) <- meta$sample_id
meta$condition <- factor(meta$condition)
if (nlevels(meta$condition) != 2) stop("Exactly two condition levels required.")
dds <- DESeqDataSetFromMatrix(countData = mat, colData = meta, design = ~ condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)
lev <- levels(meta$condition)
res <- results(dds, contrast = c("condition", lev[2], lev[1]), alpha = 0.05)
res_df <- as.data.frame(res) %>% rownames_to_column("transcript_id") %>% arrange(padj)
write_csv(res_df, file.path(out_dir, "transcript_differential_abundance.csv"))
norm <- as.data.frame(counts(dds, normalized = TRUE)) %>% rownames_to_column("transcript_id")
write_csv(norm, file.path(out_dir, "transcript_normalized_counts.csv"))
message("Transcript-feature analysis completed: ", out_dir)
message("Important: use a transcript-usage/switching method and transcript-to-gene annotation to test isoform switching.")
