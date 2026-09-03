suppressPackageStartupMessages({
  library(tidyverse)
  library(KEGGREST)
})

res <- read_tsv("DESeq2_KEGG_pathways/DESeq2_KEGG_pathways.tsv")

# pathway column = map00010
pathways <- unique(res$pathway)

get_kegg_category <- function(pid) {
  info <- tryCatch(
    keggGet(sub("map", "", pid))[[1]],
    error = function(e) NULL
  )
  
  if (is.null(info) || is.null(info$CLASS)) {
    return(tibble(
      pathway = pid,
      category_L1 = NA_character_,
      category_L2 = NA_character_
    ))
  }
  
  cls <- strsplit(info$CLASS, ";")[[1]]
  
  tibble(
    pathway = pid,
    category_L1 = trimws(cls[1]),
    category_L2 = trimws(cls[2])
  )
}

cat_map <- map_dfr(pathways, get_kegg_category)


res_cat <- res %>%
  left_join(cat_map, by = "pathway")

write_tsv(res_cat, "DESeq2_KEGG_pathways/DESeq2_KEGG_pathways_categorized.tsv")
