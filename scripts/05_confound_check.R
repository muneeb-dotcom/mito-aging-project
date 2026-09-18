library(recount3)
library(tidyverse)

tissues <- c("MUSCLE", "SKIN", "BRAIN", "LIVER", "HEART")

get_qc <- function(tissue) {
  cat("Getting QC metadata:", tissue, "\n")
  rse <- create_rse_manual(
    project = tissue,
    project_home = "data_sources/gtex",
    organism = "human",
    annotation = "gencode_v26",
    type = "gene"
  )
  tibble(
    sample_id = colnames(rse),
    tissue = tissue,
    rin = colData(rse)$gtex.smrin,
    ischemic_time = colData(rse)$gtex.smtsisch
  )
}

qc_data <- map_dfr(tissues, get_qc)
write_csv(qc_data, "data/raw/gtex/gtex_qc_metadata.csv")
cat("Saved", nrow(qc_data), "rows\n")