
suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(pheatmap)
})

card_dir    <- "/.../15_CADR/"     # directory containing CARD results (subdirectories per sample)
ko_matrix   <- "/.../09b_COUNTS_KO_matrix.tsv"  # <- PROVIDE: KO matrix (rows=KO, columns=samples) or vice versa
ko_id_col   <- "KO" # name of the KO column if the file is in long format / contains KO as a column
ko_to_path  <- "/.../11_COUNTS_KEGG_Pathway_matrix.tsv"  # <- optional: KO-to-pathway mapping table (KO, pathway_id, pathway_name)

outdir <- "/.../16_integration_KO_CADR/"
dir.create(file.path(outdir, "figures"), showWarnings = FALSE)
dir.create(file.path(outdir, "tables"), showWarnings = FALSE)

cutoff_keep <- c("Strict","Loose")

card_files <- list.files(card_dir, pattern="^card_amr_.*\\.txt$", recursive=TRUE, full.names=TRUE)
stopifnot(length(card_files) > 0)

card_all <- card_files %>%
  set_names() %>%
  map_dfr(~{
    sample_name <- basename(dirname(.x))
    read_tsv(.x, col_types = cols(.default="c"), show_col_types = FALSE) %>%
      clean_names() %>%
      mutate(sample = sample_name)
  })

must_cols <- c("sample","cut_off","aro","aro_name","best_hit_aro_category")
stopifnot(all(must_cols %in% colnames(card_all)))

card_filt <- card_all %>% filter(cut_off %in% cutoff_keep)

write_tsv(card_filt, file.path(outdir, "tables", "CARD_all_filtered.tsv"))

mech_long <- card_filt %>%
  separate_rows(best_hit_aro_category, sep=";") %>%
  mutate(best_hit_aro_category = str_trim(best_hit_aro_category)) %>%
  mutate(
    category_type = case_when(
      str_detect(best_hit_aro_category, "antibiotic inactivation") ~ "inactivation",
      str_detect(best_hit_aro_category, "efflux") ~ "efflux",
      str_detect(best_hit_aro_category, "target") ~ "target",
      str_detect(best_hit_aro_category, "reduced permeability") ~ "permeability",
      TRUE ~ "other"
    )
  )

mech_summary <- mech_long %>%
  count(sample, category_type, name="n_hits")

write_tsv(mech_summary, file.path(outdir, "tables", "AMR_mechanism_types_per_sample.tsv"))

mech_mat <- mech_summary %>%
  pivot_wider(names_from = sample, values_from = n_hits, values_fill = 0) %>%
  column_to_rownames("category_type") %>%
  as.matrix()
mech_mat_num <- as.matrix(mech_mat)

mode(mech_mat_num) <- "numeric"

mech_mat_num <- mech_mat_num[complete.cases(mech_mat_num), , drop = FALSE]

mech_mat_num <- mech_mat_num[apply(mech_mat_num, 1, var) > 0, , drop = FALSE]

mech_scaled <- t(scale(t(mech_mat_num)))   
png(file.path(outdir, "figures", "AMR_mechanism_types_heatmap.png"),
    width=1600, height=1200, res=300, bg="white")
pheatmap(mech_scaled,
         cluster_rows=TRUE, cluster_cols=TRUE,
         main="AMR mechanism TYPES (CARD)",
         fontsize_row=10, fontsize_col=10)
dev.off()

drug_long <- card_filt %>%
  separate_rows(best_hit_aro_category, sep=";") %>%
  mutate(best_hit_aro_category = str_trim(best_hit_aro_category)) %>%
  filter(str_detect(best_hit_aro_category, "resistance")) %>%
  mutate(
    drug_class = best_hit_aro_category %>%
      str_remove("^determinant of\\s+") %>%
      str_remove("\\s+resistance$") %>%
      str_trim()
  )

drug_summary <- drug_long %>%
  count(sample, drug_class, name="n_hits")

write_tsv(drug_summary, file.path(outdir, "tables", "AMR_drug_classes_per_sample.tsv"))

drug_mat <- drug_summary %>%
  pivot_wider(names_from = sample, values_from = n_hits, values_fill = 0) %>%
  column_to_rownames("drug_class") %>%
  as.matrix()
drug_mat_num <- as.matrix(drug_mat)

mode(drug_mat_num) <- "numeric"

drug_mat_num <- drug_mat_num[complete.cases(drug_mat_num), , drop = FALSE]

drug_mat_num <- drug_mat_num[apply(drug_mat_num, 1, var) > 0, , drop = FALSE]

drug_scaled <- t(scale(t(drug_mat_num)))   
png(file.path(outdir, "figures", "AMR_drug_classes_heatmap.png"),
    width=1800, height=1400, res=300, bg="white")
pheatmap(drug_scaled,
         cluster_rows=TRUE, cluster_cols=TRUE,
         main="AMR drug classes inferred from ARO categories (CARD)",
         fontsize_row=9, fontsize_col=10)
dev.off()

aro_summary <- card_filt %>%
  count(sample, aro_name, name="n_hits")

aro_mat <- aro_summary %>%
  pivot_wider(names_from=sample, values_from=n_hits, values_fill=0) %>%
  column_to_rownames("aro_name") %>%
  as.matrix()

write_tsv(aro_summary, file.path(outdir, "tables", "AMR_ARO_hits_per_sample.tsv"))

ko <- read_tsv(ko_matrix, col_types = cols(.default="c"), show_col_types = FALSE) %>% clean_names()

first_col <- colnames(ko)[1]

if (first_col %in% c("ko","k0","ko_id","koid") || str_detect(tolower(first_col), "ko")) {
  ko_wide <- ko
  colnames(ko_wide)[1] <- "ko"
} else {
  ko_wide <- ko
  colnames(ko_wide)[1] <- "ko"
}

ko_wide <- ko_wide %>%
  mutate(ko = str_remove(ko, "^ko:"))

ko_num <- ko_wide %>%
  mutate(across(-ko, ~as.numeric(.x)))  # NAs są OK

sample_cols <- setdiff(colnames(ko_num), "ko")

colnames(ko_num)[colnames(ko_num) %in% sample_cols] <- toupper(sample_cols)

card_samples <- sort(unique(card_filt$sample))
ko_samples   <- sort(sample_cols)

common_samples <- intersect(card_samples, ko_samples)
stopifnot(length(common_samples) >= 2)

message("✅ Common samples: ", paste(common_samples, collapse=", "))

drug_mat2 <- drug_mat[, common_samples, drop=FALSE]
mech_mat2 <- mech_mat[, common_samples, drop=FALSE]
aro_mat2  <- aro_mat[,  common_samples, drop=FALSE]

ko_mat <- ko_num %>%
  select(ko, all_of(common_samples)) %>%
  column_to_rownames("ko") %>%
  as.matrix()


amr_metrics <- tibble(sample = common_samples) %>%
  mutate(
    amr_total_aro_hits   = colSums(aro_mat2,  na.rm=TRUE),
    amr_total_drug_hits  = colSums(drug_mat2, na.rm=TRUE),
    amr_total_mech_hits  = colSums(mech_mat2, na.rm=TRUE),
    amr_n_unique_aro     = colSums(aro_mat2 > 0,  na.rm=TRUE),
    amr_n_unique_drugs   = colSums(drug_mat2 > 0, na.rm=TRUE),
    amr_n_unique_mech    = colSums(mech_mat2 > 0, na.rm=TRUE)
  )

write_tsv(amr_metrics, file.path(outdir, "tables", "AMR_metrics_per_sample.tsv"))

ko_metrics <- tibble(sample = common_samples) %>%
  mutate(
    ko_total_abund = colSums(ko_mat[, common_samples, drop=FALSE], na.rm=TRUE),
    ko_n_present   = colSums(ko_mat[, common_samples, drop=FALSE] > 0, na.rm=TRUE)
  )

write_tsv(ko_metrics, file.path(outdir, "tables", "KO_metrics_per_sample.tsv"))

integrated_metrics <- amr_metrics %>%
  left_join(ko_metrics, by="sample")

write_tsv(integrated_metrics, file.path(outdir, "tables", "CARD_KO_integrated_metrics.tsv"))

cor_df <- integrated_metrics %>%
  select(-sample) %>%
  cor(method="spearman", use="pairwise.complete.obs") %>%
  as.data.frame() %>%
  rownames_to_column("feature_1") %>%
  pivot_longer(-feature_1, names_to="feature_2", values_to="spearman_r")

write_tsv(cor_df, file.path(outdir, "tables", "Spearman_correlations_AMR_KO_metrics.tsv"))

cor_mat <- integrated_metrics %>%
  column_to_rownames("sample") %>%
  as.matrix()

feat_cor <- cor(cor_mat, method="spearman", use="pairwise.complete.obs")

png(file.path(outdir, "figures", "AMR_KO_metrics_correlation_heatmap.png"),
    width=1600, height=1400, res=300, bg="white")
pheatmap(feat_cor,
         cluster_rows=TRUE, cluster_cols=TRUE,
         main="Spearman correlation: AMR metrics vs KO metrics")
dev.off()
