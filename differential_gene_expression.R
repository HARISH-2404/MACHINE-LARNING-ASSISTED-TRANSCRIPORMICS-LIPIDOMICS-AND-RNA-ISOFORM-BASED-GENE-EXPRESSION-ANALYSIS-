#!/usr/bin/env Rscript
# Differential gene expression from raw integer counts using DESeq2.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript R/differential_gene_expression.R gene_counts.csv metadata.csv output_dir")
}
count_file <- args[1]; meta_file <- args[2]; out_dir <- args[3]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages({
  library(DESeq2)
  library(readr)
  library(dplyr)
  library(tibble)
})
counts_df <- read_csv(count_file, show_col_types = FALSE)
meta <- read_csv(meta_file, show_col_types = FALSE)
if (!"feature_id" %in% names(counts_df)) stop("Count matrix first column must be feature_id.")
if (!all(c("sample_id", "condition") %in% names(meta))) stop("Metadata needs sample_id and condition.")
if (anyDuplicated(counts_df$feature_id)) stop("feature_id values must be unique.")
if (anyDuplicated(meta$sample_id)) stop("sample_id values must be unique.")
sample_cols <- setdiff(names(counts_df), "feature_id")
if (!setequal(sample_cols, meta$sample_id)) stop("Sample IDs in counts and metadata do not match.")
counts_df <- counts_df %>% select(feature_id, all_of(meta$sample_id))
mat <- as.matrix(counts_df[, -1])
rownames(mat) <- counts_df$feature_id
storage.mode(mat) <- "numeric"
if (anyNA(mat) || any(mat < 0) || any(abs(mat - round(mat)) > 1e-8)) {
  stop("DESeq2 input must be non-negative, non-missing integer raw counts.")
}
storage.mode(mat) <- "integer"
meta <- as.data.frame(meta)
rownames(meta) <- meta$sample_id
meta$condition <- factor(meta$condition)
if (nlevels(meta$condition) != 2) stop("This script expects exactly two condition groups.")
if ("batch" %in% names(meta) && length(unique(na.omit(meta$batch))) > 1 && !anyNA(meta$batch)) {
  meta$batch <- factor(meta$batch)
  design_formula <- ~ batch + condition
} else {
  design_formula <- ~ condition
}
dds <- DESeqDataSetFromMatrix(countData = mat, colData = meta, design = design_formula)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)
levels_condition <- levels(meta$condition)
if (length(levels_condition) != 2) stop("Exactly two condition levels are required.")
res <- results(dds, contrast = c("condition", levels_condition[2], levels_condition[1]), alpha = 0.05)
res_df <- as.data.frame(res) %>% rownames_to_column("gene_id") %>% arrange(padj)
write_csv(res_df, file.path(out_dir, "DEG_results.csv"))
norm <- counts(dds, normalized = TRUE) %>% as.data.frame() %>% rownames_to_column("gene_id")
write_csv(norm, file.path(out_dir, "normalized_counts.csv"))
vsd <- vst(dds, blind = FALSE)
png(file.path(out_dir, "PCA.png"), width = 1200, height = 900, res = 150)
print(plotPCA(vsd, intgroup = "condition"))
dev.off()
png(file.path(out_dir, "MA_plot.png"), width = 1200, height = 900, res = 150)
plotMA(res, ylim = c(-5, 5), main = "DESeq2 MA plot")
dev.off()
message("DEG analysis completed: ", out_dir)
