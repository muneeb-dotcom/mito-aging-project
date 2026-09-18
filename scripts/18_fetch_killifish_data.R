library(GEOquery)
library(tidyverse)

# Full sample metadata across both platform subsets
gse <- getGEO("GSE66712", GSEMatrix = TRUE)
meta1 <- pData(gse[[1]])[, c("title", "geo_accession")]
meta2 <- pData(gse[[2]])[, c("title", "geo_accession")]
meta <- bind_rows(meta1, meta2)
cat("Total samples in metadata:", nrow(meta), "\n")

# Parse tissue/age/group directly from the title string
meta <- meta %>%
  mutate(
    tissue = str_extract(title, "^\\w+"),
    age_weeks = as.numeric(str_extract(title, "\\d+(?=weeks)")),
    is_normal_aging = str_detect(title, "normal aging")
  )
print(table(meta$tissue, meta$is_normal_aging))

write_csv(meta, "data/raw/killifish_sample_metadata.csv")

# Download the RPKM matrix
options(timeout = 300)
download.file(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE66nnn/GSE66712/suppl/GSE66712_nfu_rpkm_185_samples.txt.gz",
  "data/raw/GSE66712_nfu_rpkm_185_samples.txt.gz"
)

rpkm <- read.table(gzfile("data/raw/GSE66712_nfu_rpkm_185_samples.txt.gz"), header = TRUE, row.names = 1, sep = "\t")
cat("RPKM matrix dims:", dim(rpkm), "\n")
cat("First few row (gene) IDs:\n")
print(head(rownames(rpkm), 10))
cat("First few column (sample) IDs:\n")
print(head(colnames(rpkm), 10))