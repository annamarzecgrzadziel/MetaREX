#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(matrixStats)
  library(pheatmap)
})

cat("\n KO-based AMR pipeline (from scratch)\n\n")

count_file <- readline("📥 Enter the path to the KO count matrix (TSV/CSV): ")
meta_file  <- readline("📥 Enter the path to metadata (TSV/CSV): ")
outdir     <- readline("📤 Enter the output directory [AMR_results]: ")

if (outdir == "") outdir <- "AMR_results"

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "r_objects"), showWarnings = FALSE, recursive = TRUE)

read_table_auto <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") {
    read.csv(path, check.names = FALSE)
  } else {
    read.table(path, sep = "\t", header = TRUE,
               check.names = FALSE, quote = "", comment.char = "")
  }
}

counts_df <- read_table_auto(count_file)
meta_df   <- read_table_auto(meta_file)

if (!("sampleID" %in% colnames(meta_df))) {
  stop(" Metadata must contain the 'sampleID' column")
}

gene_col <- 1
ko_ids <- counts_df[[gene_col]]
counts_mat <- as.matrix(counts_df[, -gene_col, drop = FALSE])
rownames(counts_mat) <- ko_ids
mode(counts_mat) <- "numeric"
counts_mat <- round(counts_mat)

samples <- colnames(counts_mat)
meta_df <- meta_df[match(samples, meta_df$sampleID), , drop = FALSE]
rownames(meta_df) <- meta_df$sampleID

keep <- rowSums(counts_mat >= 10) >= 2
counts_filt <- counts_mat[keep, , drop = FALSE]

cat(sprintf(" Retained  wano %d / %d KO\n",
            nrow(counts_filt), nrow(counts_mat)))

dds <- DESeqDataSetFromMatrix(
  countData = counts_filt,
  colData   = meta_df,
  design    = ~ 1
)

dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind = TRUE)

saveRDS(dds, file.path(outdir, "r_objects", "dds.rds"))
saveRDS(vsd, file.path(outdir, "r_objects", "vsd.rds"))

ko_mat <- assay(vsd)

ko_amr <- c(
  # Efflux
  "K18138", "K18139", "K18140",
  "K08204", "K11744",

  # Enzymatic resistance
  "K00817", "K01467", "K03883",

  # Target modification
  "K18220", "K02542",

  # Regulators
  "K07690", "K07787", "K03466"
)

present_ko_amr <- intersect(ko_amr, rownames(ko_mat))

print(present_ko_amr)

if (length(present_ko_amr) == 0) {
  stop(" No AMR-associated KO were detected – this is also a valid result, but the pipeline ends here.")
}

amr_expr <- ko_mat[present_ko_amr, , drop = FALSE]

amr_summary <- data.frame(
  KO = rownames(amr_expr),
  mean_vst = rowMeans(amr_expr),
  sd_vst   = apply(amr_expr, 1, sd),
  stringsAsFactors = FALSE
)

amr_summary <- amr_summary[
  order(amr_summary$mean_vst, decreasing = TRUE),
]

write.table(
  amr_summary,
  file = file.path(outdir, "tables", "AMR_KO_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

amr_activity_index <- colSums(amr_expr)

write.table(
  data.frame(
    sample = names(amr_activity_index),
    AMR_activity = amr_activity_index
  ),
  file = file.path(outdir, "tables", "AMR_activity_index.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

mat_scaled <- scale(t(amr_expr))
mat_scaled[is.na(mat_scaled)] <- 0

pdf(file.path(outdir, "figures", "AMR_heatmap.pdf"), width = 7, height = 6)
pheatmap(
  mat_scaled,
  main = "Expression of antibiotic resistance–associated functions (KO-based)",
  show_rownames = TRUE,
  cluster_rows = TRUE,
  cluster_cols = TRUE
)
dev.off()

png(
  file.path(outdir, "figures", "AMR_heatmap.png"),
  width = 7,
  height = 6,
  units = "in",
  res = 300
)
pheatmap(
  mat_scaled,
  main = "Expression of antibiotic resistance–associated functions (KO-based)",
  show_rownames = TRUE,
  cluster_rows = TRUE,
  cluster_cols = TRUE
)
dev.off()

cat("\n✅ TRACK A COMPLETED SUCCESSFULLY\n")
cat("📁 Results in: ", normalizePath(outdir), "\n\n")
