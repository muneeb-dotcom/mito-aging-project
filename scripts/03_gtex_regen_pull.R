library(recount3)
library(recount)
library(tidyverse)

ecm_genes <- read_csv(r"(C:\Users\CC\regen-convergence-v2\results\de_tables\branch_b_ecm_signaling_genes.csv)", show_col_types = FALSE)$gene
core_genes <- toupper(read_csv(r"(C:\Users\CC\regen-convergence-v2\results\de_tables\convergent_genes_alias_corrected.csv)", show_col_types = FALSE)$x)
regen_genes <- unique(c(ecm_genes, core_genes))
cat("Total regen signature genes (human symbols):", length(regen_genes), "\n")

tissues <- c("MUSCLE", "SKIN", "BRAIN", "LIVER", "HEART")

pull_one <- function(tissue) {
  cat("Pulling:", tissue, "\n")
  rse <- create_rse_manual(
    project = tissue,
    project_home = "data_sources/gtex",
    organism = "human",
    annotation = "gencode_v26",
    type = "gene"
  )
  assay(rse, "counts") <- transform_counts(rse)
  assays(rse)$TPM <- getTPM(rse, length_var = "bp_length")

  gene_symbols <- rowData(rse)$gene_name
  keep <- gene_symbols %in% regen_genes
  rse <- rse[keep, ]
  gene_symbols <- gene_symbols[keep]

  tpm <- as.data.frame(assay(rse, "TPM"))
  tpm$gene <- gene_symbols

  tpm %>%
    pivot_longer(-gene, names_to = "sample_id", values_to = "tpm") %>%
    mutate(tissue = tissue)
}

all_data <- map_dfr(tissues, pull_one)
write_csv(all_data, "data/raw/gtex/gtex_regen_signature_long.csv")
cat("Saved", nrow(all_data), "rows,", n_distinct(all_data$gene), "unique genes matched\n")