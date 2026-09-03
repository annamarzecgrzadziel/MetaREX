#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ALDEx2)
  library(tidyverse)
  library(ggplot2)
})

counts_file <- "/.../11_COUNTS_KEGG_Pathway_matrix.tsv"
group_file  <- readline(" Enter the sample groups file (TSV/CSV: sampleID, group): ")
out_dir     <- "/.../10a_DEG_KEGG_Pathway_ALDEx2"

pairwise_dir <- file.path(out_dir, "pairwise_tables")
plots_dir    <- file.path(out_dir, "plots")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(pairwise_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

mc_samples   <- 128
fdr_cutoff   <- 0.05
effect_cutoff <- 1

cat("   mc.samples:", mc_samples, "\n")
cat("   FDR cutoff:", fdr_cutoff, "\n")
cat("   effect cutoff:", effect_cutoff, "\n\n")

read_table_auto <- function(path) {
  if (!file.exists(path)) {
    stop(" File not found: ", path)
  }
  
  first_line <- readLines(path, n = 1, warn = FALSE)
  
  if (grepl("\t", first_line)) {
    readr::read_tsv(path, show_col_types = FALSE)
  } else if (grepl(";", first_line)) {
    readr::read_delim(path, delim = ";", show_col_types = FALSE)
  } else {
    readr::read_csv(path, show_col_types = FALSE)
  }
}

is_valid_tsv <- function(path, required_cols = character(), min_rows = 1) {
  if (!file.exists(path) || file.info(path)$size <= 0) {
    return(FALSE)
  }
  
  tab <- tryCatch(
    readr::read_tsv(path, show_col_types = FALSE, progress = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(tab) || nrow(tab) < min_rows) {
    return(FALSE)
  }
  
  all(required_cols %in% colnames(tab))
}

save_plot_if_missing <- function(plot, pdf_file, png_file, width, height) {
  if (!file.exists(pdf_file) || file.info(pdf_file)$size <= 0) {
    ggsave(
      filename = pdf_file,
      plot = plot,
      width = width,
      height = height
    )
    cat("    PDF:", basename(pdf_file), "\n")
  } else {
    cat("    PDF already exists:", basename(pdf_file), "\n")
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
    cat("    PNG:", basename(png_file), "\n")
  } else {
    cat("    PNG already exists:", basename(png_file), "\n")
  }
}

run_aldex_pair <- function(count_mat, group_vec, pair_name) {
  cat("\n ALDEx2:", pair_name, "\n")
  print(table(group_vec))
  
  if (length(unique(group_vec)) != 2) {
    stop(" Comparison ", pair_name, " does not contain exactly two groups.")
  }
  
  if (ncol(count_mat) != length(group_vec)) {
    stop(" The number of samples does not match the group vector length for ", pair_name)
  }
  
  keep_features <- rowSums(count_mat) > 0
  count_mat <- count_mat[keep_features, , drop = FALSE]
  
  cat("   Pathways after removing zero-count features:", nrow(count_mat), "\n")
  
  if (nrow(count_mat) == 0) {
    stop(" No KEGG Pathways remain after removing zero-count features.")
  }
  
  clr <- aldex.clr(
    reads = count_mat,
    conds = group_vec,
    mc.samples = mc_samples,
    denom = "all",
    verbose = FALSE
  )
  
  tt <- as.data.frame(
    aldex.ttest(clr, paired.test = FALSE)
  )
  
  eff <- as.data.frame(
    aldex.effect(clr)
  )
  
  tt$feature <- rownames(tt)
  eff$feature <- rownames(eff)
  
  tt %>%
    select(feature, we.ep, we.eBH, wi.ep, wi.eBH) %>%
    left_join(
      eff %>%
        select(feature, diff.btw, diff.win, effect, overlap),
      by = "feature"
    ) %>%
    mutate(
      comparison = pair_name,
      group_1 = sort(unique(group_vec))[1],
      group_2 = sort(unique(group_vec))[2],
      signif = we.eBH < fdr_cutoff & abs(effect) > effect_cutoff,
      status = if_else(signif, "significant", "not_significant"),
      neglog10_padj = -log10(we.eBH + 1e-300)
    ) %>%
    select(
      feature,
      comparison,
      group_1,
      group_2,
      we.ep,
      we.eBH,
      wi.ep,
      wi.eBH,
      diff.btw,
      diff.win,
      effect,
      overlap,
      signif,
      status,
      neglog10_padj
    )
}

counts <- read_table_auto(counts_file)

if (ncol(counts) < 2) {
  stop(" The COUNTS matrix must contain an ID column and at least one sample.")
}

feature_col <- colnames(counts)[1]

if (anyDuplicated(counts[[feature_col]]) > 0) {
  duplicated_ids <- unique(counts[[feature_col]][duplicated(counts[[feature_col]])])
  stop(
    " Duplicate KEGG Pathway identifiers: ",
    paste(head(duplicated_ids, 20), collapse = ", ")
  )
}

count_mat <- counts %>%
  column_to_rownames(feature_col) %>%
  as.data.frame()

count_mat[] <- lapply(
  count_mat,
  function(x) suppressWarnings(as.numeric(as.character(x)))
)

count_mat <- as.matrix(count_mat)

if (anyNA(count_mat)) {
  stop(" NA values appeared after converting COUNTS to numeric.")
}

if (any(count_mat < 0)) {
  stop(" The COUNTS matrix contains negative values.")
}

non_integer_n <- sum(count_mat %% 1 != 0, na.rm = TRUE)

if (non_integer_n > 0) {
  stop(
    " The matrix contains ", non_integer_n,
    " non-integer values. ALDEx2 requires raw counts."
  )
}

count_mat <- round(count_mat)
storage.mode(count_mat) <- "integer"

cat("", nrow(count_mat), "KEGG Pathways ×", ncol(count_mat), "próbek\n")

groups_df <- read_table_auto(group_file)

required_meta_cols <- c("sampleID", "group")
missing_meta_cols <- setdiff(required_meta_cols, colnames(groups_df))

if (length(missing_meta_cols) > 0) {
  stop(
    " Missing metadata columns: ",
    paste(missing_meta_cols, collapse = ", "),
    ". Required: sampleID, group."
  )
}

groups_df <- groups_df %>%
  transmute(
    sampleID = as.character(sampleID),
    group = as.character(group)
  )

if (anyNA(groups_df$sampleID) || any(groups_df$sampleID == "")) {
  stop(" Metadata contains an empty sampleID.")
}

if (anyNA(groups_df$group) || any(groups_df$group == "")) {
  stop(" Metadata contains an empty group value.")
}

if (anyDuplicated(groups_df$sampleID) > 0) {
  stop(" Metadata contains duplicate sampleID values.")
}

missing_in_counts <- setdiff(groups_df$sampleID, colnames(count_mat))
missing_in_meta   <- setdiff(colnames(count_mat), groups_df$sampleID)

if (length(missing_in_counts) > 0) {
  stop(
    " Samples from metadata not found in the COUNTS matrix: ",
    paste(missing_in_counts, collapse = ", ")
  )
}

if (length(missing_in_meta) > 0) {
  cat(
    " The COUNTS matrix contains samples without metadata; they will be excluded:\n   ",
    paste(missing_in_meta, collapse = ", "),
    "\n"
  )
}

count_mat <- count_mat[, groups_df$sampleID, drop = FALSE]
groups <- groups_df$group
names(groups) <- groups_df$sampleID

if (ncol(count_mat) != length(groups)) {
  stop(" The number of COUNTS columns does not match the number of groups.")
}

print(table(groups))

unique_groups <- unique(groups)

if (length(unique_groups) < 2) {
  stop(" At least two groups are required for comparative analysis.")
}

group_pairs <- combn(unique_groups, 2, simplify = FALSE)

required_pair_cols <- c(
  "feature", "comparison", "group_1", "group_2",
  "we.ep", "we.eBH", "wi.ep", "wi.eBH",
  "diff.btw", "diff.win", "effect", "overlap",
  "signif", "status", "neglog10_padj"
)

aldex_pairwise_results <- list()

for (pair in group_pairs) {
  pair_name <- paste(pair, collapse = "_vs_")
  pair_file <- file.path(
    pairwise_dir,
    paste0("ALDEx2_KEGG_Pathway_", pair_name, ".tsv")
  )
  
  if (is_valid_tsv(pair_file, required_pair_cols, min_rows = 1)) {
    cat("\n Skipping completed comparison:", pair_name, "\n")
    aldex_pairwise_results[[pair_name]] <- read_tsv(
      pair_file,
      show_col_types = FALSE,
      progress = FALSE
    )
    next
  }
  
  if (file.exists(pair_file)) {
    cat(" Existing file is empty or incomplete — recalculating:", pair_name, "\n")
  }
  
  keep_samples <- groups %in% pair
  count_mat_sub <- count_mat[, keep_samples, drop = FALSE]
  
  groups_sub <- factor(groups[keep_samples], levels = pair)
  groups_sub <- as.character(groups_sub)
  
  result <- run_aldex_pair(
    count_mat = count_mat_sub,
    group_vec = groups_sub,
    pair_name = pair_name
  )
  
  write_tsv(result, pair_file)
  cat(" Saved comparison:", pair_file, "\n")
  
  aldex_pairwise_results[[pair_name]] <- result
}

aldex_pairwise_df <- bind_rows(aldex_pairwise_results)

if (nrow(aldex_pairwise_df) == 0) {
  stop(" No comparison results were obtained.")
}

all_results_file <- file.path(
  out_dir,
  "ALDEx2_pairwise_all_comparisons.tsv"
)

write_tsv(aldex_pairwise_df, all_results_file)
cat("\n Saved/refreshed combined results table:\n", all_results_file, "\n")

cat("\n Volcano plot for all comparisons...\n")

p_all <- ggplot(
  aldex_pairwise_df,
  aes(x = effect, y = neglog10_padj, color = status)
) +
  geom_point(alpha = 0.7, size = 2) +
  geom_hline(
    yintercept = -log10(fdr_cutoff),
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_vline(
    xintercept = c(-effect_cutoff, effect_cutoff),
    linetype = "dashed",
    linewidth = 0.4
  ) +
  facet_wrap(~ comparison, scales = "free") +
  theme_bw(base_size = 12) +
  labs(
    title = "ALDEx2 KEGG Pathway volcano plots",
    x = "Effect size",
    y = "-log10(adjusted p-value)",
    color = "Status"
  )

save_plot_if_missing(
  plot = p_all,
  pdf_file = file.path(plots_dir, "ALDEx2_volcano_all_comparisons.pdf"),
  png_file = file.path(plots_dir, "ALDEx2_volcano_all_comparisons.png"),
  width = 14,
  height = 10
)

cat("\n Individual volcano plots...\n")

for (cmp in unique(aldex_pairwise_df$comparison)) {
  df_sub <- aldex_pairwise_df %>%
    filter(comparison == cmp)
  
  p <- ggplot(
    df_sub,
    aes(x = effect, y = neglog10_padj, color = status)
  ) +
    geom_point(alpha = 0.7, size = 2) +
    geom_hline(
      yintercept = -log10(fdr_cutoff),
      linetype = "dashed",
      linewidth = 0.4
    ) +
    geom_vline(
      xintercept = c(-effect_cutoff, effect_cutoff),
      linetype = "dashed",
      linewidth = 0.4
    ) +
    theme_bw(base_size = 12) +
    labs(
      title = paste("ALDEx2 KEGG Pathway:", cmp),
      x = "Effect size",
      y = "-log10(adjusted p-value)",
      color = "Status"
    )
  
  save_plot_if_missing(
    plot = p,
    pdf_file = file.path(
      plots_dir,
      paste0("ALDEx2_volcano_", cmp, ".pdf")
    ),
    png_file = file.path(
      plots_dir,
      paste0("ALDEx2_volcano_", cmp, ".png")
    ),
    width = 8,
    height = 6
  )
}

cat("\n Summary of significant pathway counts...\n")

summary_df <- aldex_pairwise_df %>%
  group_by(comparison) %>%
  summarise(
    n_total = n(),
    n_significant = sum(signif, na.rm = TRUE),
    n_positive_effect = sum(signif & effect > 0, na.rm = TRUE),
    n_negative_effect = sum(signif & effect < 0, na.rm = TRUE),
    .groups = "drop"
  )

summary_file <- file.path(
  out_dir,
  "ALDEx2_significant_feature_counts.tsv"
)

write_tsv(summary_df, summary_file)
cat(" Saved summary:", summary_file, "\n")

p_sum <- ggplot(
  summary_df,
  aes(x = comparison, y = n_significant)
) +
  geom_col() +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "Number of significant KEGG Pathways per comparison",
    x = "Comparison",
    y = "Number of significant KEGG Pathways"
  )

save_plot_if_missing(
  plot = p_sum,
  pdf_file = file.path(plots_dir, "ALDEx2_significant_feature_counts.pdf"),
  png_file = file.path(plots_dir, "ALDEx2_significant_feature_counts.png"),
  width = 8,
  height = 5
)

if (length(unique_groups) > 2) {
  cat("\n Global ALDEx2 Kruskal-Wallis test...\n")
  
  global_file <- file.path(
    out_dir,
    "ALDEx2_KEGG_Pathway_global_KW.tsv"
  )
  
  excluded_file <- file.path(
    out_dir,
    "ALDEx2_KEGG_Pathway_global_KW_excluded_features.tsv"
  )
  
  required_kw_cols <- c("feature", "kw.ep", "kw.eBH")
  
  if (is_valid_tsv(global_file, required_kw_cols, min_rows = 1)) {
    cat(" A valid global KW test result already exists — skipping.\n")
  } else {
    keep_global <- rowSums(count_mat) > 0
    global_mat <- count_mat[keep_global, , drop = FALSE]
    
    cat("   Input pathways:", nrow(count_mat), "\n")
    cat("   Pathways after removing zero-count features:", nrow(global_mat), "\n")
    
    aldex_global <- aldex.clr(
      reads = global_mat,
      conds = groups,
      mc.samples = mc_samples,
      denom = "all",
      verbose = FALSE
    )
    
    kw <- as.data.frame(aldex.kw(aldex_global))
    kw$feature <- rownames(kw)
    
    kw_raw <- kw %>%
      select(feature, any_of(c("kw.ep", "kw.eBH")))
    
    if (!all(required_kw_cols %in% colnames(kw_raw))) {
      stop(
        " aldex.kw() did not return the expected columns: kw.ep and kw.eBH."
      )
    }
    
    kw_res <- tibble(
      feature = rownames(count_mat)
    ) %>%
      left_join(kw_raw, by = "feature") %>%
      arrange(kw.eBH)
    
    write_tsv(kw_res, global_file)
    cat(" Saved global KW test:\n", global_file, "\n")
    
    excluded <- kw_res %>%
      filter(is.na(kw.ep) | is.na(kw.eBH))
    
    if (nrow(excluded) > 0) {
      write_tsv(excluded, excluded_file)
      cat(
        " ALDEx2 excluded ", nrow(excluded),
        " pathways. List:\n", excluded_file, "\n",
        sep = ""
      )
    } else if (file.exists(excluded_file)) {
      file.remove(excluded_file)
    }
  }
}


cat("\n Results saved in:\n")
cat(out_dir, "\n\n")

cat("Main output files:\n")
cat(" - ALDEx2_pairwise_all_comparisons.tsv\n")
cat(" - ALDEx2_significant_feature_counts.tsv\n")
cat(" - pairwise_tables/\n")
cat(" - plots/\n")

if (length(unique_groups) > 2) {
  cat(" - ALDEx2_KEGG_Pathway_global_KW.tsv\n")
}

cat("\n STEP 10 KEGG PATHWAY COMPLETED\n")
