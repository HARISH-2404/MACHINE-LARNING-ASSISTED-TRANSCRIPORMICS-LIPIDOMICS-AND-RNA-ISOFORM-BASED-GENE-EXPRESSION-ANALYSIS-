#!/usr/bin/env Rscript
# Import Salmon quant.sf files and test gene-level differential expression with tximport + DESeq2.
# Usage: Rscript R/tximport_deseq2.R samples.csv tx2gene.csv quant_root output_dir
# samples.csv columns: sample_id,condition,batch(optional)
# tx2gene.csv columns: transcript_id,gene_id; IDs must match the indexed transcriptome.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4) stop("Usage: Rscript R/tximport_deseq2.R samples.csv tx2gene.csv quant_root output_dir")
suppressPackageStartupMessages({library(tximport); library(DESeq2); library(readr); library(dplyr); library(tibble)})
meta <- read_csv(args[1], show_col_types=FALSE)
tx2gene <- read_csv(args[2], show_col_types=FALSE)
quant_root <- args[3]; out <- args[4]; dir.create(out, recursive=TRUE, showWarnings=FALSE)
stopifnot(all(c("sample_id","condition") %in% names(meta)))
stopifnot(all(c("transcript_id","gene_id") %in% names(tx2gene)))
if (anyDuplicated(meta$sample_id)) stop("sample_id must be unique")
files <- file.path(quant_root, meta$sample_id, "quant.sf")
if (!all(file.exists(files))) stop("Missing quant.sf files: ", paste(files[!file.exists(files)], collapse=", "))
names(files) <- meta$sample_id
meta <- as.data.frame(meta); rownames(meta) <- meta$sample_id; meta$condition <- factor(meta$condition)
if (nlevels(meta$condition) != 2) stop("This script expects exactly two conditions")
txi <- tximport(files, type="salmon", tx2gene=as.data.frame(tx2gene), ignoreTxVersion=TRUE)
if ("batch" %in% names(meta) && !anyNA(meta$batch) && length(unique(meta$batch)) > 1) {
  meta$batch <- factor(meta$batch); design <- ~ batch + condition
} else design <- ~ condition
dds <- DESeqDataSetFromTximport(txi, colData=meta, design=design)
dds <- dds[rowSums(counts(dds)) >= 10, ]; dds <- DESeq(dds)
lev <- levels(meta$condition)
res <- results(dds, contrast=c("condition", lev[2], lev[1]), alpha=0.05)
out_df <- as.data.frame(res) %>% rownames_to_column("gene_id") %>% arrange(padj)
write_csv(out_df, file.path(out,"DEG_tximport_DESeq2.csv"))
write_csv(as.data.frame(counts(dds, normalized=TRUE)) %>% rownames_to_column("gene_id"), file.path(out,"normalized_gene_counts.csv"))
png(file.path(out,"MA_plot.png"), width=1200,height=900,res=150); plotMA(res, ylim=c(-5,5)); dev.off()
message("tximport/DESeq2 complete. Confirm reference and transcript annotation compatibility.")