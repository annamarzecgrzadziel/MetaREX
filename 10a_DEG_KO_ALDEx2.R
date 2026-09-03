#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ALDEx2)
  library(tidyverse)
  library(ggplot2)
})


ko_file   <- readline("📥 Enter the KO COUNTS matrix file, e.g. 09b_COUNTS_KO_matrix.tsv: ")
meta_file <- readline("📥 Enter the metadata file (sampleID → group): ")
out_dir   <- readline("📤 Enter the output directory: ")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "pairwise_tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "plots"), showWarnings = FALSE, recursive = TRUE)

mc_samples <- 128
fdr_cutoff <- 0.05
effect_cutoff <- 1

cat("   mc.samples:", mc_samples, "\n")
cat("   FDR cutoff:", fdr_cutoff, "\n")
cat("   effect cutoff:", effect_cutoff, "\n\n")

is_nonempty_file <- function(path) {
  file.exists(path) && !is.na(file.info(path)$size) && file.info(path)$size > 0
}

is_valid_tsv <- function(path, required_cols = NULL) {
  if (!is_nonempty_file(path)) {
    return(FALSE)
  }
  
  ok <- tryCatch({
    x <- readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
    
    if (nrow(x) == 0) {
      return(FALSE)
    }
    
    if (!is.null(required_cols) && !all(required_cols %in% colnames(x))) {
      return(FALSE)
    }
    
    TRUE
  }, error = function(e) {
    FALSE
  })
  
  isTRUE(ok)
}

safe_ggsave <- function(filename, plot, width, height, dpi = NULL) {
  if (is_nonempty_file(filename)) {
    cat("  Valid plot file already exists – skipping:", filename, "\n")
    return(invisible(FALSE))
  }
  
  args <- list(
    filename = filename,
    plot = plot,
    width = width,
    height = height
  )
  
  if (!is.null(dpi)) {
    args$dpi <- dpi
  }
  
  do.call(ggsave, args)
  invisible(TRUE)
}

pair_required_cols <- c(
  "KO", "comparison", "group_1", "group_2",
  "we.ep", "we.eBH", "wi.ep", "wi.eBH",
  "diff.btw", "diff.win", "effect", "overlap",
  "neglog10_FDR", "change", "significant"
)

ko <- read_tsv(ko_file, show_col_types = FALSE)
meta <- read_tsv(meta_file, show_col_types = FALSE)

required_meta_cols <- c("sampleID", "group")

missing_meta_cols <- setdiff(required_meta_cols, colnames(meta))

if (length(missing_meta_cols) > 0) {
  stop(
    " Missing metadata columns: ",
    paste(missing_meta_cols, collapse = ", "),
    "\nMetadata must contain the columns: sampleID and group"
  )
}

meta <- meta %>%
  mutate(
    sampleID = as.character(sampleID),
    group = as.character(group)
  )

if (any(is.na(meta$sampleID)) || any(meta$sampleID == "")) {
  stop(" The sampleID column contains empty values.")
}

if (any(is.na(meta$group)) || any(meta$group == "")) {
  stop(" The group column contains empty values.")
}

ko_id_col <- colnames(ko)[1]

cat(" First column in the KO matrix treated as feature ID:", ko_id_col, "\n")

missing_samples <- setdiff(meta$sampleID, colnames(ko))

if (length(missing_samples) > 0) {
  stop(
    " The following samples from metadata are missing from the KO matrix:\n",
    paste(missing_samples, collapse = ", ")
  )
}

ko_mat <- ko %>%
  column_to_rownames(ko_id_col) %>%
  select(all_of(meta$sampleID)) %>%
  as.data.frame()

ko_mat[] <- lapply(ko_mat, function(x) as.numeric(as.character(x)))
ko_mat <- as.matrix(ko_mat)

if (any(is.na(ko_mat))) {
  stop(" NA values appeared after conversion to numeric. Check whether the KO matrix contains text in sample columns.")
}

cat("   mode:", mode(ko_mat), "\n")
cat("   value range:\n")
print(summary(as.vector(ko_mat)))

non_integer_n <- sum(ko_mat %% 1 != 0, na.rm = TRUE)

cat("   number of non-integer values:", non_integer_n, "\n")

if (non_integer_n > 0) {
  cat("\n The matrix contains non-integer values.\n")
  cat("   Examples of non-integer values:\n")
  print(head(as.vector(ko_mat)[ko_mat %% 1 != 0], 20))
   stop(
    "\nALDEx2 requires raw counts as integer values.\n",
    "You probably provided a TPM file, e.g. from 09a_aggregate_tpm_by_eggnog.\n",
    "Use a counts file instead, e.g. 09b_COUNTS_KO_matrix.tsv.\n"
  )
}

ko_mat <- round(ko_mat)
storage.mode(ko_mat) <- "integer"

groups <- meta$group

if (ncol(ko_mat) != length(groups)) {
  stop(" The number of columns in the KO matrix does not match the length of the groups vector.")
}

cat("\n Loaded:", nrow(ko_mat), "KO ×", ncol(ko_mat), "samples\n")
cat(" Groups:\n")
print(table(groups))

run_aldex_pair <- function(count_mat, group_vec, pair_name = "comparison") {
  
  cat("\n ALDEx2:", pair_name, "\n")
  print(table(group_vec))
  
  if (length(unique(group_vec)) != 2) {
    stop(" run_aldex_pair requires exactly two groups.")
  }
  
  keep_features <- rowSums(count_mat) > 0
  count_mat <- count_mat[keep_features, , drop = FALSE]
  
  cat("   KO after removing zero-count features:", nrow(count_mat), "\n")
  
  aldex <- aldex.clr(
    reads = count_mat,
    conds = group_vec,
    mc.samples = mc_samples,
    denom = "all",
    verbose = FALSE
  )
  
  tt <- aldex.ttest(
    aldex,
    paired.test = FALSE
  )
  
  effect <- aldex.effect(aldex)
  
  tt <- as.data.frame(tt)
  effect <- as.data.frame(effect)
  
  tt$KO <- rownames(tt)
  effect$KO <- rownames(effect)
  
  res <- tt %>%
    select(KO, we.ep, we.eBH, wi.ep, wi.eBH) %>%
    left_join(
      effect %>%
        select(KO, diff.btw, diff.win, effect, overlap),
      by = "KO"
    ) %>%
    mutate(
      comparison = pair_name,
      group_1 = sort(unique(group_vec))[1],
      group_2 = sort(unique(group_vec))[2],
      neglog10_FDR = -log10(we.eBH + 1e-300),
      change = case_when(
        we.eBH < fdr_cutoff & effect >= effect_cutoff  ~ paste0("Higher_in_", group_2),
        we.eBH < fdr_cutoff & effect <= -effect_cutoff ~ paste0("Higher_in_", group_1),
        TRUE ~ "Stable"
      ),
      significant = we.eBH < fdr_cutoff & abs(effect) >= effect_cutoff
    ) %>%
    select(
      KO,
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
      neglog10_FDR,
      change,
      significant
    )
  
  return(res)
}

unique_groups <- unique(groups)
aldex_results <- list()

if (length(unique_groups) == 2) {
  group_pairs <- list(unique_groups)
  cat(" Exactly 2 groups detected – one comparison.\n")
} else {
  group_pairs <- combn(unique_groups, 2, simplify = FALSE)
  cat(" More than 2 groups detected – pairwise comparisons.\n")
}

for (pair in group_pairs) {
  
  pair_name <- paste(pair, collapse = "_vs_")
  
  out_pair_tsv <- file.path(
    out_dir,
    "pairwise_tables",
    paste0("10_ALDEx2_KO_", pair_name, ".tsv")
  )
  
  if (is_valid_tsv(out_pair_tsv, pair_required_cols)) {
    cat("\n  Comparison result already exists – loading and skipping calculations:", pair_name, "\n")
    aldex_results[[pair_name]] <- read_tsv(
      out_pair_tsv,
      show_col_types = FALSE,
      progress = FALSE
    )
    next
  }
  
  if (file.exists(out_pair_tsv)) {
    cat("\n Existing file is empty, corrupted, or incomplete – recalculating:", pair_name, "\n")
  }
  
  keep_samples <- groups %in% pair
  count_mat_sub <- ko_mat[, keep_samples, drop = FALSE]
  groups_sub <- factor(groups[keep_samples], levels = pair)
  groups_sub <- as.character(groups_sub)
  
  pair_result <- run_aldex_pair(
    count_mat = count_mat_sub,
    group_vec = groups_sub,
    pair_name = pair_name
  )
  
  write_tsv(pair_result, out_pair_tsv)
  cat(" Saved comparison result:", out_pair_tsv, "\n")
  
  aldex_results[[pair_name]] <- pair_result
}

aldex_df <- bind_rows(aldex_results)
out_all_tsv <- file.path(out_dir, "10_ALDEx2_KO_all_comparisons.tsv")

write_tsv(aldex_df, out_all_tsv)

cat(" Saved/rebuilt combined table:\n", out_all_tsv, "\n")
cat(" Comparison tables are located in:\n", file.path(out_dir, "pairwise_tables"), "\n")

cat("\n Generating volcano plots...\n")

for (cmp in unique(aldex_df$comparison)) {
  
  df_sub <- aldex_df %>%
    filter(comparison == cmp)
  
  p <- ggplot(df_sub, aes(x = effect, y = neglog10_FDR, color = change)) +
    geom_point(alpha = 0.75, size = 2) +
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
      title = paste("ALDEx2 KO volcano plot:", cmp),
      x = "Effect size",
      y = "-log10(FDR)",
      color = "Change"
    )
  
  ggsave(
    filename = file.path(out_dir, "plots", paste0("10_ALDEx2_KO_volcano_", cmp, ".pdf")),
    plot = p,
    width = 8,
    height = 6
  )
  
  ggsave(
    filename = file.path(out_dir, "plots", paste0("10_ALDEx2_KO_volcano_", cmp, ".png")),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
}

cat(" Generating TOP effect size plots...\n")

for (cmp in unique(aldex_df$comparison)) {
  
  top <- aldex_df %>%
    filter(comparison == cmp) %>%
    filter(significant) %>%
    arrange(desc(abs(effect))) %>%
    slice_head(n = 30)
  
  if (nrow(top) == 0) {
    cat(" No significant KO for comparison:", cmp, "- skipping TOP effect plot.\n")
    next
  }
  
  p <- ggplot(top, aes(x = reorder(KO, effect), y = effect, fill = change)) +
    geom_col() +
    coord_flip() +
    theme_bw(base_size = 12) +
    labs(
      title = paste("Top differential KO:", cmp),
      x = "KO",
      y = "Effect size",
      fill = "Change"
    )
  
  ggsave(
    filename = file.path(out_dir, "plots", paste0("10_ALDEx2_KO_top_effect_", cmp, ".pdf")),
    plot = p,
    width = 8,
    height = 10
  )
  
  ggsave(
    filename = file.path(out_dir, "plots", paste0("10_ALDEx2_KO_top_effect_", cmp, ".png")),
    plot = p,
    width = 8,
    height = 10,
    dpi = 300
  )
}

cat(" Generating summary of significant KO counts...\n")

summary_df <- aldex_df %>%
  group_by(comparison) %>%
  summarise(
    n_total = n(),
    n_significant = sum(significant, na.rm = TRUE),
    n_positive_effect = sum(significant & effect > 0, na.rm = TRUE),
    n_negative_effect = sum(significant & effect < 0, na.rm = TRUE),
    .groups = "drop"
  )

write_tsv(
  summary_df,
  file.path(out_dir, "10_ALDEx2_KO_summary_counts.tsv")
)

p_sum <- ggplot(summary_df, aes(x = comparison, y = n_significant)) +
  geom_col() +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "Number of significant KO per comparison",
    x = "Comparison",
    y = "Number of significant KO"
  )

safe_ggsave(
  filename = file.path(out_dir, "plots", "10_ALDEx2_KO_significant_counts.pdf"),
  plot = p_sum,
  width = 8,
  height = 5
)

safe_ggsave(
  filename = file.path(out_dir, "plots", "10_ALDEx2_KO_significant_counts.png"),
  plot = p_sum,
  width = 8,
  height = 5,
  dpi = 300
)

if (length(unique_groups) > 2) {
  
  out_kw <- file.path(out_dir, "10_ALDEx2_KO_global_KW.tsv")
  out_kw_excluded <- file.path(out_dir, "10_ALDEx2_KO_global_KW_excluded_features.tsv")
  
  if (is_valid_tsv(out_kw, c("KO", "kw.ep", "kw.eBH"))) {
    cat("\n  Valid global KW test result already exists – skipping calculations:\n", out_kw, "\n")
  } else {
    
    if (file.exists(out_kw)) {
      cat("\n Existing KW result is empty, corrupted, or incomplete – recalculating.\n")
    }
    
    cat("\n Running global ALDEx2 Kruskal-Wallis test for multiple groups...\n")
    
    global_keep <- rowSums(ko_mat) > 0
    ko_mat_global <- ko_mat[global_keep, , drop = FALSE]
    
     cat("   Input KO:", nrow(ko_mat), "\n")
    cat("   KO after removing zero-count features:", nrow(ko_mat_global), "\n")
    
    aldex_global <- aldex.clr(
      reads = ko_mat_global,
      conds = groups,
      mc.samples = mc_samples,
      denom = "all",
      verbose = FALSE
    )
    
    kw <- aldex.kw(aldex_global)
    kw_raw <- as.data.frame(kw) %>%
      tibble::rownames_to_column("KO") %>%
      select(KO, any_of(c("kw.ep", "kw.eBH")))
    
    if (!all(c("kw.ep", "kw.eBH") %in% colnames(kw_raw))) {
      stop(" aldex.kw did not return the required columns kw.ep and kw.eBH.")
    }
    
    kw_res <- tibble(KO = rownames(ko_mat)) %>%
      left_join(kw_raw, by = "KO") %>%
      arrange(is.na(kw.eBH), kw.eBH)
    
    excluded_kw <- kw_res %>%
      filter(is.na(kw.ep) | is.na(kw.eBH))
    
    cat("   KO with test result:", nrow(kw_res) - nrow(excluded_kw), "\n")
    cat("   KO excluded by the test:", nrow(excluded_kw), "\n")
    
    write_tsv(kw_res, out_kw)
    
    if (nrow(excluded_kw) > 0) {
      write_tsv(excluded_kw, out_kw_excluded)
      cat(" Lista KO pominiętych przez test została zapisana w:\n", out_kw_excluded, "\n")
    } else if (file.exists(out_kw_excluded)) {
      file.remove(out_kw_excluded)
    }
    
    cat(" Saved global KW test:\n", out_kw, "\n")
  }
}

cat("\n Output files saved in:\n")
cat(out_dir, "\n\n")

cat("Main output files:\n")
cat(" - 10_ALDEx2_KO_all_comparisons.tsv\n")
cat(" - 10_ALDEx2_KO_summary_counts.tsv\n")
cat(" - pairwise_tables/\n")
cat(" - plots/\n")
