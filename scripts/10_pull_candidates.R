library(recount3)
library(recount)
library(tidyverse)

candidate_genes <- c("NAMPT", "SIRT1", "PPARGC1A")
tissues <- c("MUSCLE", "SKIN", "BRAIN", "LIVER", "HEART")

pull_one <- function(tissue) {
  cat("Pulling:", tissue, "\n")
  rse <- create_rse_manual(
    project = tissue, project_home = "data_sources/gtex",
    organism = "human", annotation = "gencode_v26", type = "gene"
  )
  assay(rse, "counts") <- transform_counts(rse)
  assays(rse)$TPM <- getTPM(rse, length_var = "bp_length")

  gene_symbols <- rowData(rse)$gene_name
  keep <- gene_symbols %in% candidate_genes
  rse <- rse[keep, ]
  gene_symbols <- gene_symbols[keep]

  tpm <- as.data.frame(assay(rse, "TPM"))
  tpm$gene <- gene_symbols

  tpm %>%
    pivot_longer(-gene, names_to = "sample_id", values_to = "tpm") %>%
    mutate(tissue = tissue)
}

all_data <- map_dfr(tissues, pull_one)
write_csv(all_data, "data/raw/gtex/gtex_candidates_long.csv")
cat("Saved", nrow(all_data), "rows,", n_distinct(all_data$gene), "genes matched\n")