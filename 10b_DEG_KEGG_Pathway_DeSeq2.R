#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(apeglm)
  library(pheatmap)
  library(ggplot2)
})

counts_file <- "kegg_pathway_counts.tsv"
meta_file   <- "metadata.tsv"
out_dir     <- "DESeq2_KEGG_Pathway_results"

ref_level     <- "control"
min_counts    <- 20
alpha_cutoff  <- 0.05
lfc_cutoff    <- 1
top_heatmap_n <- 50

pairwise_dir <- file.path(out_dir, "pairwise_tables")
plots_dir    <- file.path(out_dir, "plots")
state_dir    <- file.path(out_dir, ".state")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(pairwise_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(state_dir, showWarnings = FALSE, recursive = TRUE)

read_table_auto <- function(path) {

  if (!file.exists(path)) {
    stop("File not found: ", path)
  }

  first_line <- readLines(path, n = 1, warn = FALSE)

  if (grepl("\t", first_line)) {

    readr::read_tsv(
      path,
      show_col_types = FALSE
    )

  } else if (grepl(";", first_line)) {

    readr::read_delim(
      path,
      delim = ";",
      show_col_types = FALSE
    )

  } else {

    readr::read_csv(
      path,
      show_col_types = FALSE
    )
  }
}


is_valid_tsv <- function(
    path,
    required_cols = character(),
    min_rows = 1
) {

  if (!file.exists(path) || file.info(path)$size <= 0) {
    return(FALSE)
  }

  tab <- tryCatch(
    readr::read_tsv(
      path,
      show_col_types = FALSE,
      progress = FALSE
    ),
    error = function(e) NULL
  )

  if (is.null(tab)) {
    return(FALSE)
  }

  if (nrow(tab) < min_rows) {
    return(FALSE)
  }

  all(required_cols %in% colnames(tab))
}


is_valid_rds <- function(path) {

  if (!file.exists(path) || file.info(path)$size <= 0) {
    return(FALSE)
  }

  obj <- tryCatch(
    readRDS(path),
    error = function(e) NULL
  )

  !is.null(obj)
}


safe_name <- function(x) {

  gsub(
    "[^A-Za-z0-9_.-]+",
    "_",
    x
  )
}


save_plot_if_missing <- function(
    plot,
    pdf_file,
    png_file,
    width,
    height
) {

  if (!file.exists(pdf_file) || file.info(pdf_file)$size <= 0) {

    ggsave(
      filename = pdf_file,
      plot = plot,
      width = width,
      height = height
    )

    cat("   PDF:", basename(pdf_file), "\n")

  } else {

    cat("   PDF already exists:", basename(pdf_file), "\n")
  }


  if (!file.exists(png_file) || file.info(png_file)$size <= 0) {

    ggsave(
      filename = png_file,
      plot = plot,
      width = width,
      height = height,
      dpi = 300,
      device = "png",
      type = "cairo"
    )

    cat("   PNG:", basename(png_file), "\n")

  } else {

    cat("   PNG already exists:", basename(png_file), "\n")
  }
}


build_input_signature <- function(
    counts_file,
    meta_file,
    design_formula,
    ref_level,
    min_counts
) {

  paste(
    normalizePath(counts_file),
    file.info(counts_file)$size,
    as.character(file.info(counts_file)$mtime),
    normalizePath(meta_file),
    file.info(meta_file)$size,
    as.character(file.info(meta_file)$mtime),
    deparse(design_formula),
    ref_level,
    min_counts,
    sep = "|"
  )
}

counts <- read_table_auto(counts_file)

if (!"KEGG_Pathway" %in% colnames(counts)) {

  stop(
    "The counts matrix must contain a 'KEGG_Pathway' column."
  )
}

if (anyDuplicated(counts$KEGG_Pathway) > 0) {

  duplicated_ids <- unique(
    counts$KEGG_Pathway[duplicated(counts$KEGG_Pathway)]
  )

  stop(
    "Duplicate KEGG_Pathway identifiers: ",
    paste(
      head(duplicated_ids, 20),
      collapse = ", "
    )
  )
}

counts_mat <- counts %>%
  column_to_rownames("KEGG_Pathway") %>%
  as.data.frame()

counts_mat[] <- lapply(
  counts_mat,
  function(x) {
    suppressWarnings(
      as.numeric(as.character(x))
    )
  }
)

counts_mat <- as.matrix(counts_mat)

if (anyNA(counts_mat)) {

  stop(
    " NA values appeared after converting counts to numeric."
  )
}

if (any(counts_mat < 0)) {

  stop(
    " The counts matrix contains negative values."
  )
}

non_integer_n <- sum(
  counts_mat %% 1 != 0,
  na.rm = TRUE
)

if (non_integer_n > 0) {

  stop(
    " DESeq2 requires raw counts. ",
    "Number of non-integer values: ",
    non_integer_n
  )
}

counts_mat <- round(counts_mat)
storage.mode(counts_mat) <- "integer"

cat(
  " ",
  nrow(counts_mat),
  " KEGG Pathway × ",
  ncol(counts_mat),
  " próbek\n",
  sep = ""
)

meta <- read_table_auto(meta_file)

required_meta_cols <- c(
  "sample",
  "condition"
)

missing_meta_cols <- setdiff(
  required_meta_cols,
  colnames(meta)
)

if (length(missing_meta_cols) > 0) {

  stop(
    " Missing metadata columns: ",
    paste(
      missing_meta_cols,
      collapse = ", "
    ),
    ". Required columns: sample, condition."
  )
}

meta <- meta %>%
  mutate(
    sample = as.character(sample),
    condition = as.character(condition)
  )

if (anyNA(meta$sample) || any(meta$sample == "")) {

  stop(
    " Metadata contains an empty sample value."
  )
}

if (anyNA(meta$condition) || any(meta$condition == "")) {

  stop(
    " Metadata contains an empty condition value."
  )
}

if (anyDuplicated(meta$sample) > 0) {

  duplicated_samples <- unique(
    meta$sample[duplicated(meta$sample)]
  )

  stop(
    " Duplicate samples in metadata: ",
    paste(
      head(duplicated_samples, 20),
      collapse = ", "
    )
  )
}

missing_in_meta <- setdiff(
  colnames(counts_mat),
  meta$sample
)

missing_in_counts <- setdiff(
  meta$sample,
  colnames(counts_mat)
)

if (length(missing_in_meta) > 0) {

  stop(
    " Missing metadata for samples from the counts matrix: ",
    paste(
      missing_in_meta,
      collapse = ", "
    )
  )
}

if (length(missing_in_counts) > 0) {

  cat(
    " Samples present in metadata but absent from counts will be omitted:\n   ",
    paste(
      missing_in_counts,
      collapse = ", "
    ),
    "\n"
  )
}

meta <- meta %>%
  filter(sample %in% colnames(counts_mat)) %>%
  column_to_rownames("sample")

meta <- meta[
  colnames(counts_mat),
  ,
  drop = FALSE
]

if (!identical(
  rownames(meta),
  colnames(counts_mat)
)) {

  stop(
    " Failed to correctly match sample order."
  )
}

meta$condition <- factor(meta$condition)

if (!ref_level %in% levels(meta$condition)) {

  stop(
    " Reference level '",
    ref_level,
    "' nie występuje w condition.\n",
    "Dostępne poziomy: ",
    paste(
      levels(meta$condition),
      collapse = ", "
    )
  )
}

meta$condition <- relevel(
  meta$condition,
  ref = ref_level
)

has_batch <- "batch" %in% colnames(meta)

if (has_batch) {

  meta$batch <- factor(meta$batch)

  if (anyNA(meta$batch)) {

    stop(
      " The batch column contains NA values."
    )
  }

  if (nlevels(meta$batch) < 2) {

    cat(
      " The batch column has only one level — ",
      "batch will be omitted from the design.\n"
    )

    has_batch <- FALSE
  }
}

design_formula <- if (has_batch) {
  ~ batch + condition
} else {
  ~ condition
}

cat(
  "✅ Design: ",
  deparse(design_formula),
  "\n",
  sep = ""
)

cat("✅ Groups:\n")
print(table(meta$condition))

if (nlevels(meta$condition) < 2) {

  stop(
    " At least two groups are required for comparative analysis."
  )
}

dds_file <- file.path(
  state_dir,
  "dds_fitted.rds"
)

vsd_file <- file.path(
  state_dir,
  "vsd.rds"
)

signature_file <- file.path(
  state_dir,
  "input_signature.txt"
)

input_signature <- build_input_signature(
  counts_file = counts_file,
  meta_file = meta_file,
  design_formula = design_formula,
  ref_level = ref_level,
  min_counts = min_counts
)

state_matches <- FALSE

if (file.exists(signature_file)) {

  saved_signature <- readLines(
    signature_file,
    warn = FALSE
  )

  state_matches <- identical(
    saved_signature,
    input_signature
  )
}

dds <- NULL

if (
  state_matches &&
  is_valid_rds(dds_file)
) {

  candidate_dds <- readRDS(dds_file)

  if (
    identical(
      colnames(candidate_dds),
      colnames(counts_mat)
    )
  ) {

    dds <- candidate_dds

    cat(
      "\n Loading previously fitted DESeq2 object.\n"
    )

  } else {

    cat(
      " Saved DESeq2 object does not match the current samples — ",
      "recalculating.\n"
    )
  }
}

if (is.null(dds)) {

  dds <- DESeqDataSetFromMatrix(
    countData = counts_mat,
    colData = meta,
    design = design_formula
  )

  keep_features <- rowSums(
    counts(dds)
  ) >= min_counts

  cat(
    " Prefiltering: retaining ",
    sum(keep_features),
    " of ",
    length(keep_features),
    " KEGG Pathways\n",
    sep = ""
  )

  dds <- dds[
    keep_features,
  ]

  if (nrow(dds) == 0) {

    stop(
      " No KEGG Pathways remain after prefiltering."
    )
  }

  dds <- DESeq(dds)

  saveRDS(
    dds,
    dds_file
  )

  writeLines(
    input_signature,
    signature_file
  )

  if (file.exists(vsd_file)) {
    file.remove(vsd_file)
  }

  cat(
    " Saved DESeq2 object: ",
    dds_file,
    "\n",
    sep = ""
  )
}

if (
  state_matches &&
  is_valid_rds(vsd_file)
) {

  cat(
    " Loading previously calculated VST.\n"
  )

  vsd <- readRDS(vsd_file)

} else {

  cat(" calculating VST...\n")

  if (nrow(dds) >= 1000) {

    vsd <- vst(
      dds,
      blind = FALSE
    )

  } else {

    cat(
      " Fewer than 1000 features (",
      nrow(dds),
      ") — using varianceStabilizingTransformation() instead of vst().\n",
      sep = ""
    )

    vsd <- varianceStabilizingTransformation(
      dds,
      blind = FALSE
    )
  }

  saveRDS(
    vsd,
    vsd_file
  )
}

pca_pdf <- file.path(
  plots_dir,
  "PCA_KEGG_Pathway_DESeq2.pdf"
)

pca_png <- file.path(
  plots_dir,
  "PCA_KEGG_Pathway_DESeq2.png"
)

if (
  !file.exists(pca_pdf) ||
  file.info(pca_pdf)$size <= 0 ||
  !file.exists(pca_png) ||
  file.info(pca_png)$size <= 0
) {

  intgroup_cols <- if (has_batch) {
    c("condition", "batch")
  } else {
    "condition"
  }

  pca <- plotPCA(
    vsd,
    intgroup = intgroup_cols
  ) +
    ggtitle("PCA – KEGG Pathway (DESeq2)") +
    theme_bw(base_size = 12)

  save_plot_if_missing(
    plot = pca,
    pdf_file = pca_pdf,
    png_file = pca_png,
    width = 7,
    height = 6
  )

} else {

  cat(
    " PCA already exists — skipping.\n"
  )
}

condition_levels <- levels(
  meta$condition
)

group_pairs <- combn(
  condition_levels,
  2,
  simplify = FALSE
)

required_result_cols <- c(
  "KEGG_Pathway",
  "comparison",
  "group_1",
  "group_2",
  "baseMean",
  "log2FoldChange",
  "lfcSE",
  "pvalue",
  "padj",
  "status",
  "significant",
  "neglog10_padj"
)

pairwise_results <- list()

for (pair in group_pairs) {

  group_1 <- pair[1]
  group_2 <- pair[2]

  comparison <- paste(
    group_1,
    group_2,
    sep = "_vs_"
  )

  comparison_safe <- safe_name(
    comparison
  )

  result_file <- file.path(
    pairwise_dir,
    paste0(
      "DESeq2_KEGG_Pathway_",
      comparison_safe,
      ".tsv"
    )
  )

  if (
    is_valid_tsv(
      result_file,
      required_result_cols,
      min_rows = 1
    )
  ) {

    cat(
      " Skipping completed comparison: ",
      comparison,
      "\n",
      sep = ""
    )

    pairwise_results[[comparison]] <- read_tsv(
      result_file,
      show_col_types = FALSE,
      progress = FALSE
    )

    next
  }

  if (file.exists(result_file)) {

    cat(
      " Existing result is empty or incomplete — recalculating: ",
      comparison,
      "\n",
      sep = ""
    )

  } else {

    cat(
      " Comparison: ",
      comparison,
      "\n",
      sep = ""
    )
  }

  contrast_vec <- c(
    "condition",
    group_2,
    group_1
  )

  if (group_1 == ref_level) {

    expected_coef <- paste0(
      "condition_",
      make.names(group_2),
      "_vs_",
      make.names(group_1)
    )

    matching_coef <- resultsNames(dds)[
      resultsNames(dds) == expected_coef
    ]

    if (length(matching_coef) == 1) {

      res_shr <- lfcShrink(
        dds,
        coef = matching_coef,
        type = "apeglm"
      )

    } else {

      cat(
        " Could not find an unambiguous coef for apeglm for ",
        comparison,
        ". Używam type='normal'.\n",
        sep = ""
      )

      res_shr <- lfcShrink(
        dds,
        contrast = contrast_vec,
        type = "normal"
      )
    }

  } else {

    res_shr <- lfcShrink(
      dds,
      contrast = contrast_vec,
      type = "normal"
    )
  }

  res_tbl <- as.data.frame(
    res_shr
  ) %>%
    rownames_to_column("KEGG_Pathway") %>%
    mutate(
      comparison = comparison,
      group_1 = group_1,
      group_2 = group_2,

      significant = (
        !is.na(padj) &
          padj < alpha_cutoff &
          abs(log2FoldChange) > lfc_cutoff
      ),

      status = case_when(

        significant &
          log2FoldChange > 0 ~
          paste0(
            "Higher_in_",
            group_2
          ),

        significant &
          log2FoldChange < 0 ~
          paste0(
            "Higher_in_",
            group_1
          ),

        TRUE ~ "NS"
      ),

      neglog10_padj = -log10(
        padj + 1e-300
      )
    ) %>%
    select(
      KEGG_Pathway,
      comparison,
      group_1,
      group_2,
      baseMean,
      log2FoldChange,
      lfcSE,
      pvalue,
      padj,
      status,
      significant,
      neglog10_padj
    ) %>%
    arrange(padj)

  write_tsv(
    res_tbl,
    result_file
  )

  cat(
    " Written: ",
    result_file,
    "\n",
    sep = ""
  )

  pairwise_results[[comparison]] <- res_tbl
}

all_results <- bind_rows(
  pairwise_results
)

if (nrow(all_results) == 0) {

  stop(
    " No DESeq2 results were obtained."
  )
}

all_results_file <- file.path(
  out_dir,
  "DESeq2_KEGG_Pathway_all_comparisons.tsv"
)

write_tsv(
  all_results,
  all_results_file
)

cat(
  "\n Saved/refreshed combined results table:\n",
  all_results_file,
  "\n"
)

cat("\n📈 MA plots...\n")

for (comparison in names(pairwise_results)) {

  comparison_safe <- safe_name(
    comparison
  )

  ma_pdf <- file.path(
    plots_dir,
    paste0(
      "MAplot_KEGG_Pathway_DESeq2_",
      comparison_safe,
      ".pdf"
    )
  )

  ma_png <- file.path(
    plots_dir,
    paste0(
      "MAplot_KEGG_Pathway_DESeq2_",
      comparison_safe,
      ".png"
    )
  )

  if (
    file.exists(ma_pdf) &&
    file.info(ma_pdf)$size > 0 &&
    file.exists(ma_png) &&
    file.info(ma_png)$size > 0
  ) {

    cat(
      " MA plot already exists: ",
      comparison,
      "\n",
      sep = ""
    )

    next
  }

  ma_data <- pairwise_results[[comparison]] %>%
    filter(
      !is.na(baseMean),
      !is.na(log2FoldChange)
    )

  p_ma <- ggplot(
    ma_data,
    aes(
      x = baseMean,
      y = log2FoldChange,
      color = significant
    )
  ) +
    geom_point(
      alpha = 0.6,
      size = 1.5
    ) +
    scale_x_log10() +
    geom_hline(
      yintercept = 0,
      linetype = "dashed"
    ) +
    coord_cartesian(
      ylim = c(-6, 6)
    ) +
    theme_bw(base_size = 12) +
    labs(
      title = paste(
        "MA plot – KEGG Pathway:",
        comparison
      ),
      x = "Mean normalized count",
      y = "Shrunken log2 fold change",
      color = "Significant"
    )

  save_plot_if_missing(
    plot = p_ma,
    pdf_file = ma_pdf,
    png_file = ma_png,
    width = 7,
    height = 6
  )
}

cat("\n🌋 Volcano plots...\n")

for (comparison in names(pairwise_results)) {

  comparison_safe <- safe_name(
    comparison
  )

  volcano_pdf <- file.path(
    plots_dir,
    paste0(
      "Volcano_KEGG_Pathway_DESeq2_",
      comparison_safe,
      ".pdf"
    )
  )

  volcano_png <- file.path(
    plots_dir,
    paste0(
      "Volcano_KEGG_Pathway_DESeq2_",
      comparison_safe,
      ".png"
    )
  )

  if (
    file.exists(volcano_pdf) &&
    file.info(volcano_pdf)$size > 0 &&
    file.exists(volcano_png) &&
    file.info(volcano_png)$size > 0
  ) {

    cat(
      " Volcano plot already exists: ",
      comparison,
      "\n",
      sep = ""
    )

    next
  }

  res_tbl <- pairwise_results[[comparison]]

  volcano <- ggplot(
    res_tbl,
    aes(
      x = log2FoldChange,
      y = neglog10_padj,
      color = status
    )
  ) +
    geom_point(
      alpha = 0.7,
      size = 2
    ) +
    geom_vline(
      xintercept = c(
        -lfc_cutoff,
        lfc_cutoff
      ),
      linetype = "dashed"
    ) +
    geom_hline(
      yintercept = -log10(alpha_cutoff),
      linetype = "dashed"
    ) +
    theme_bw(base_size = 12) +
    labs(
      title = paste(
        "Volcano – KEGG Pathway (DESeq2):",
        comparison
      ),
      x = paste0(
        "Shrunken log2FC (",
        res_tbl$group_2[1],
        " vs ",
        res_tbl$group_1[1],
        ")"
      ),
      y = "-log10(adjusted p-value)",
      color = "Status"
    )

  save_plot_if_missing(
    plot = volcano,
    pdf_file = volcano_pdf,
    png_file = volcano_png,
    width = 7,
    height = 6
  )
}

all_volcano_pdf <- file.path(
  plots_dir,
  "Volcano_KEGG_Pathway_DESeq2_all_comparisons.pdf"
)

all_volcano_png <- file.path(
  plots_dir,
  "Volcano_KEGG_Pathway_DESeq2_all_comparisons.png"
)

if (
  !file.exists(all_volcano_pdf) ||
  file.info(all_volcano_pdf)$size <= 0 ||
  !file.exists(all_volcano_png) ||
  file.info(all_volcano_png)$size <= 0
) {

  p_all_volcano <- ggplot(
    all_results,
    aes(
      x = log2FoldChange,
      y = neglog10_padj,
      color = status
    )
  ) +
    geom_point(
      alpha = 0.65,
      size = 1.5
    ) +
    geom_vline(
      xintercept = c(
        -lfc_cutoff,
        lfc_cutoff
      ),
      linetype = "dashed"
    ) +
    geom_hline(
      yintercept = -log10(alpha_cutoff),
      linetype = "dashed"
    ) +
    facet_wrap(
      ~ comparison,
      scales = "free"
    ) +
    theme_bw(base_size = 11) +
    labs(
      title = "DESeq2 KEGG Pathway volcano plots",
      x = "Shrunken log2 fold change",
      y = "-log10(adjusted p-value)",
      color = "Status"
    )

  save_plot_if_missing(
    plot = p_all_volcano,
    pdf_file = all_volcano_pdf,
    png_file = all_volcano_png,
    width = 14,
    height = 10
  )

} else {

  cat(
    " Combined volcano plot already exists — skipping.\n"
  )
}

heatmap_pdf <- file.path(
  plots_dir,
  "Heatmap_significant_KEGG_Pathway_DESeq2.pdf"
)

heatmap_png <- file.path(
  plots_dir,
  "Heatmap_significant_KEGG_Pathway_DESeq2.png"
)

significant_pathways <- all_results %>%
  filter(significant) %>%
  arrange(padj) %>%
  distinct(KEGG_Pathway) %>%
  slice_head(
    n = top_heatmap_n
  ) %>%
  pull(KEGG_Pathway)

if (length(significant_pathways) > 1) {

  if (
    !file.exists(heatmap_pdf) ||
    file.info(heatmap_pdf)$size <= 0 ||
    !file.exists(heatmap_png) ||
    file.info(heatmap_png)$size <= 0
  ) {

    significant_pathways <- intersect(
      significant_pathways,
      rownames(vsd)
    )

    heat_mat <- assay(vsd)[
      significant_pathways,
      ,
      drop = FALSE
    ]

    annotation_col <- as.data.frame(
      colData(vsd)
    ) %>%
      select(
        any_of(
          c(
            "condition",
            "batch"
          )
        )
      )

    pdf(
      heatmap_pdf,
      width = 10,
      height = 10
    )

    pheatmap(
      heat_mat,
      scale = "row",
      annotation_col = annotation_col,
      show_rownames = TRUE,
      fontsize_row = 7,
      main = "Top significant KEGG Pathway – DESeq2"
    )

    dev.off()

    png(
      heatmap_png,
      width = 3000,
      height = 3000,
      res = 300,
      type = "cairo"
    )

    pheatmap(
      heat_mat,
      scale = "row",
      annotation_col = annotation_col,
      show_rownames = TRUE,
      fontsize_row = 7,
      main = "Top significant KEGG Pathway – DESeq2"
    )

    dev.off()

    cat(
      " Saved heatmap of significant KEGG Pathways.\n"
    )

  } else {

    cat(
      " Heatmap already exists — skipping.\n"
    )
  }

} else {

  cat(
    " Too few significant KEGG Pathways to generate a heatmap.\n"
  )
}

summary_df <- all_results %>%
  group_by(
    comparison,
    group_1,
    group_2
  ) %>%
  summarise(
    n_total = n(),

    n_tested = sum(
      !is.na(padj)
    ),

    n_significant = sum(
      significant,
      na.rm = TRUE
    ),

    n_higher_group_1 = sum(
      significant &
        log2FoldChange < 0,
      na.rm = TRUE
    ),

    n_higher_group_2 = sum(
      significant &
        log2FoldChange > 0,
      na.rm = TRUE
    ),

    .groups = "drop"
  )

summary_file <- file.path(
  out_dir,
  "DESeq2_KEGG_Pathway_summary_counts.tsv"
)

write_tsv(
  summary_df,
  summary_file
)

p_summary <- ggplot(
  summary_df,
  aes(
    x = comparison,
    y = n_significant
  )
) +
  geom_col() +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  ) +
  labs(
    title = "Number of significant KEGG Pathway per comparison",
    x = "Comparison",
    y = "Number of significant KEGG Pathway"
  )

save_plot_if_missing(
  plot = p_summary,
  pdf_file = file.path(
    plots_dir,
    "DESeq2_KEGG_Pathway_summary_counts.pdf"
  ),
  png_file = file.path(
    plots_dir,
    "DESeq2_KEGG_Pathway_summary_counts.png"
  ),
  width = 8,
  height = 5
)

cat("\n DESeq2 KEGG Pathway analysis finished\n")
cat(" Results: ", out_dir, "\n", sep = "")
cat(" Combined table: DESeq2_KEGG_Pathway_all_comparisons.tsv\n")
cat(" Summary: DESeq2_KEGG_Pathway_summary_counts.tsv\n")
cat(" Pairwise results: pairwise_tables/\n")
cat(" Plots: plots/\n")
cat(" Resume objects: .state/\n\n")
