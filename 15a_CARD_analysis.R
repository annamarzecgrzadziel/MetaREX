
library(tidyverse)
library(janitor)
library(pheatmap)
library(FactoMineR)
library(factoextra)

card_dir <- "15_CARD"     # parent directory containing CARD results
cutoff_to_keep <- c("Strict", "Loose")

outdir <- "16_CARD_analysis"
dir.create(outdir, showWarnings = FALSE)
dir.create(file.path(outdir, "figures"), showWarnings = FALSE)
dir.create(file.path(outdir, "tables"), showWarnings = FALSE)

card_files <- list.files(
  card_dir,
  pattern = "^card_amr_.*\\.txt$",
  recursive = TRUE,
  full.names = TRUE
)

stopifnot(length(card_files) > 0)

message("🔍 Found CARD files: ", length(card_files))

card_all <- card_files %>%
  set_names() %>%
  map_dfr(~ {
    sample_name <- basename(dirname(.x))
    
    read_tsv(.x, show_col_types = FALSE) %>%
      clean_names() %>%
      mutate(sample = sample_name)
  })

card_filt <- card_all %>%
  filter(cut_off %in% cutoff_to_keep)

write_tsv(card_filt,
          file.path(outdir, "tables", "CARD_all_filtered.tsv"))

aro_summary <- card_filt %>%
  count(sample, aro_name, sort = TRUE)

write_tsv(aro_summary,
          file.path(outdir, "tables", "CARD_ARO_per_sample.tsv"))

mech_summary <- card_filt %>%
  separate_rows(aro_category, sep = ";") %>%
  filter(str_detect(aro_category, "resistance")) %>%
  count(sample, aro_category)

write_tsv(mech_summary,
          file.path(outdir, "tables", "CARD_mechanisms_per_sample.tsv"))

mech_mat <- mech_summary %>%
  pivot_wider(
    names_from = sample,
    values_from = n,
    values_fill = 0
  ) %>%
  column_to_rownames("aro_category") %>%
  as.matrix()

write_tsv(as.data.frame(mech_mat) %>%
            rownames_to_column("mechanism"),
          file.path(outdir, "tables", "CARD_mechanism_matrix.tsv"))

pheatmap(
  scale(mech_mat),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "AMR mechanisms across samples",
  fontsize_row = 8,
  fontsize_col = 10,
  filename = file.path(outdir, "figures", "AMR_mechanisms_heatmap.pdf")
)

pca_res <- PCA(
  t(mech_mat),
  scale.unit = TRUE,
  graph = FALSE
)

pdf(file.path(outdir, "figures", "AMR_PCA.pdf"), width = 7, height = 6)
fviz_pca_ind(
  pca_res,
  geom = "point",
  repel = TRUE,
  title = "PCA of AMR profiles (CARD)"
)
dev.off()

drug_class_summary <- card_filt %>%
  separate_rows(drug_class, sep = ";") %>%
  count(sample, drug_class, sort = TRUE)

write_tsv(drug_class_summary,
          file.path(outdir, "tables", "CARD_drug_classes_per_sample.tsv"))

