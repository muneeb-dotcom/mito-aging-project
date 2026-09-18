library(recount3)
library(recount)
library(tidyverse)

nuclear_oxphos <- read_csv("data/raw/mitocarta/nuclear_oxphos_genes.csv", show_col_types = FALSE)$gene_symbol
mtdna_genes <- c("MT-ND1","MT-ND2","MT-ND3","MT-ND4","MT-ND4L","MT-ND5","MT-ND6",
                  "MT-CO1","MT-CO2","MT-CO3","MT-CYB","MT-ATP6","MT-ATP8")
target_genes <- c(nuclear_oxphos, mtdna_genes)

tissues <- c("MUSCLE", "SKIN", "BRAIN", "LIVER", "HEART")

age_to_mid <- function(bracket) {
  parts <- str_split_fixed(bracket, "-", 2)
  (as.numeric(parts[,1]) + as.numeric(parts[,2])) / 2
}

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
  keep <- gene_symbols %in% target_genes
  rse <- rse[keep, ]
  gene_symbols <- gene_symbols[keep]

  tpm <- as.data.frame(assay(rse, "TPM"))
  tpm$gene <- gene_symbols

  tpm %>%
    pivot_longer(-gene, names_to = "sample_id", values_to = "tpm") %>%
    mutate(
      tissue = tissue,
      age_bracket = colData(rse)$gtex.age[match(sample_id, colnames(rse))],
      age_mid = age_to_mid(age_bracket)
    )
}

all_data <- map_dfr(tissues, pull_one)
write_csv(all_data, "data/raw/gtex/gtex_oxphos_long.csv")
cat("Saved", nrow(all_data), "rows to data/raw/gtex/gtex_oxphos_long.csv\n")