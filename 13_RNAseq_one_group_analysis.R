#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(matrixStats)
  library(pheatmap)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: 13_RNAseq_one_group_analysis_KO.R COUNTS_TSV METADATA_TSV OUTPUT_DIR")
}

count_file <- args[1]
meta_file <- args[2]
outdir <- args[3]

dir.create(file.path(outdir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "r_objects"), recursive = TRUE, showWarnings = FALSE)

read_tsv_base <- function(path) {
  read.table(path, sep = "\t", header = TRUE, check.names = FALSE,
             quote = "", comment.char = "", stringsAsFactors = FALSE)
}

scale_rows <- function(x) {
  x_scaled <- t(scale(t(x)))
  x_scaled[is.na(x_scaled)] <- 0
  x_scaled
}

counts_df <- read_tsv_base(count_file)
meta_df <- read_tsv_base(meta_file)

if (ncol(counts_df) < 3) stop("The count matrix must contain one feature column and at least two samples.")
if (!"sample" %in% colnames(meta_df)) {
  stop("Metadata must contain a 'sample' column matching count-matrix columns.")
}

feature_col <- colnames(counts_df)[1]
feature_ids <- as.character(counts_df[[feature_col]])
if (anyDuplicated(feature_ids)) stop("The feature identifiers in the count matrix are not unique.")

counts_mat <- as.matrix(counts_df[, -1, drop = FALSE])
mode(counts_mat) <- "numeric"
if (anyNA(counts_mat) || any(counts_mat < 0) || any(counts_mat %% 1 != 0)) {
  stop("The KO matrix must contain non-negative integer counts without NA values.")
}
rownames(counts_mat) <- feature_ids

if (anyDuplicated(meta_df$sample)) stop("Sample identifiers in metadata are not unique.")
if (!setequal(colnames(counts_mat), meta_df$sample)) {
  stop("Sample identifiers do not match between the count matrix and metadata.")
}
meta_df <- meta_df[match(colnames(counts_mat), meta_df$sample), , drop = FALSE]
rownames(meta_df) <- meta_df$sample

keep <- rowSums(counts_mat >= 10) >= 2
counts_filt <- counts_mat[keep, , drop = FALSE]
if (nrow(counts_filt) < 2) stop("Fewer than two KO features remain after filtering.")

dds <- DESeqDataSetFromMatrix(countData = counts_filt, colData = meta_df, design = ~ 1)
dds <- estimateSizeFactors(dds, type = "poscounts")
vsd <- vst(dds, blind = TRUE)
mat <- assay(vsd)

saveRDS(dds, file.path(outdir, "r_objects", "dds_KO.rds"))
saveRDS(vsd, file.path(outdir, "r_objects", "vsd_KO.rds"))

pca <- prcomp(t(mat), scale. = FALSE)
percent_var <- 100 * pca$sdev^2 / sum(pca$sdev^2)
pca_df <- data.frame(sample = rownames(pca$x), PC1 = pca$x[, 1], PC2 = pca$x[, 2])
if ("group" %in% colnames(meta_df)) {
  pca_df$group <- meta_df[pca_df$sample, "group"]
}

p_pca <- ggplot(pca_df, aes(PC1, PC2, label = sample)) +
  geom_point(size = 4) + geom_text(vjust = -0.7) + theme_bw() +
  xlab(sprintf("PC1: %.1f%%", percent_var[1])) +
  ylab(sprintf("PC2: %.1f%%", percent_var[2]))
ggsave(file.path(outdir, "figures", "PCA_vst_KO.pdf"), p_pca, width = 7, height = 6)
ggsave(file.path(outdir, "figures", "PCA_vst_KO.png"), p_pca, width = 7, height = 6, dpi = 300)
write.table(pca_df, file.path(outdir, "tables", "PCA_coordinates_KO.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

row_vars <- rowVars(mat)
make_heatmap <- function(n_top) {
  n_top <- min(n_top, nrow(mat))
  idx <- order(row_vars, decreasing = TRUE)[seq_len(n_top)]
  mat_top <- scale_rows(mat[idx, , drop = FALSE])
  pheatmap(mat_top, show_rownames = TRUE, cluster_rows = FALSE,
           cluster_cols = TRUE, fontsize_row = 6, fontsize_col = 10,
           main = sprintf("KO VST: top %d features by variance", n_top),
           filename = file.path(outdir, "figures", sprintf("Heatmap_vst_KO_top%d.pdf", n_top)),
           width = 8, height = max(8, n_top * 0.12))
}
make_heatmap(50)
make_heatmap(100)

summary_df <- data.frame(
  KO = rownames(mat), mean_vst = rowMeans(mat), sd_vst = rowSds(mat),
  variance_vst = row_vars
)
summary_df <- summary_df[order(summary_df$mean_vst, decreasing = TRUE), ]
write.table(summary_df, file.path(outdir, "tables", "KO_summary_mean_sd_variance.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(head(summary_df, 100), file.path(outdir, "tables", "Top100_KO_by_mean_vst.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

writeLines(c(
  paste("count_file", normalizePath(count_file), sep = "\t"),
  paste("metadata_file", normalizePath(meta_file), sep = "\t"),
  paste("n_input_KO", nrow(counts_mat), sep = "\t"),
  paste("n_retained_KO", nrow(counts_filt), sep = "\t"),
  "filter\tcount >= 10 in at least 2 samples",
  "transformation\tDESeq2 vst(blind=TRUE)",
  "heatmap_row_clustering\tFALSE",
  "heatmap_column_clustering\tTRUE"
), file.path(outdir, "analysis_parameters.tsv"))

message("KO-based VST/PCA/heatmap analysis completed: ", normalizePath(outdir))
