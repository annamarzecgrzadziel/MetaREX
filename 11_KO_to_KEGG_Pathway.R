#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: 11_KO_to_KEGG_Pathway_deduplicated.R KO_COUNTS_TSV EGGNOG_ANNOTATIONS OUTPUT_DIR")
}

ko_file <- args[1]
eggnog_file <- args[2]
out_dir <- args[3]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

normalize_pathway <- function(x) {
  x <- trimws(x)
  x <- sub("^(ko|map)", "", x, ignore.case = TRUE)
  paste0("map", x)
}

ko_mat <- read_tsv(ko_file, show_col_types = FALSE)
if (ncol(ko_mat) < 3) stop("The KO count matrix must contain a KO column and at least two sample columns.")
ko_mat <- rename(ko_mat, KO = 1)
if (anyDuplicated(ko_mat$KO)) stop("The KO count matrix contains duplicated KO identifiers.")
sample_cols <- colnames(ko_mat)[-1]

annotation_lines <- readLines(eggnog_file, warn = FALSE)
header_line <- annotation_lines[grepl("^#query", annotation_lines)][1]
if (is.na(header_line) || !length(header_line)) stop("The eggNOG annotation file has no #query header.")
header <- strsplit(sub("^#", "", header_line), "\t", fixed = FALSE)[[1]]

eggnog <- read_tsv(eggnog_file, comment = "#", col_names = header,
                   col_types = cols(.default = "c"), show_col_types = FALSE)
required <- c("KEGG_ko", "KEGG_Pathway")
if (!all(required %in% colnames(eggnog))) stop("The eggNOG file lacks KEGG_ko or KEGG_Pathway.")

ko_pathway <- eggnog %>%
  select(KEGG_ko, KEGG_Pathway) %>%
  filter(!is.na(KEGG_ko), KEGG_ko != "", !is.na(KEGG_Pathway), KEGG_Pathway != "") %>%
  separate_rows(KEGG_ko, sep = ",") %>%
  separate_rows(KEGG_Pathway, sep = ",") %>%
  mutate(
    KEGG_ko = sub("^ko:", "", trimws(KEGG_ko), ignore.case = TRUE),
    pathway_id = normalize_pathway(KEGG_Pathway)
  ) %>%
  filter(grepl("^K[0-9]{5}$", KEGG_ko), grepl("^map[0-9]+$", pathway_id)) %>%
  distinct(KEGG_ko, pathway_id)

ko_long <- ko_mat %>%
  pivot_longer(all_of(sample_cols), names_to = "sample", values_to = "count")

pathway_long <- ko_long %>%
  inner_join(ko_pathway, by = c("KO" = "KEGG_ko")) %>%
  group_by(pathway_id, sample) %>%
  summarise(count = sum(as.numeric(count)), .groups = "drop")

pathway_mat <- pathway_long %>%
  pivot_wider(names_from = sample, values_from = count, values_fill = 0) %>%
  arrange(pathway_id) %>%
  rename(KEGG_Pathway = pathway_id)

write_tsv(pathway_mat, file.path(out_dir, "11_COUNTS_KEGG_Pathway_matrix_deduplicated.tsv"))
write_tsv(ko_pathway, file.path(out_dir, "KO_to_Pathway_mapping_normalized.tsv"))
write_tsv(tibble(
  input_ko_rows = nrow(ko_mat),
  input_pathway_rows = n_distinct(ko_pathway$pathway_id),
  output_pathway_rows = nrow(pathway_mat),
  ko_pathway_mapping_rows = nrow(ko_pathway)
), file.path(out_dir, "pathway_deduplication_summary.tsv"))

message("Created ", nrow(pathway_mat), " unique KEGG Pathway features.")
