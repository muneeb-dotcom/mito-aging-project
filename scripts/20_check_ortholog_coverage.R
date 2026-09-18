library(tidyverse)

rpkm <- read.table(gzfile("data/raw/GSE66712_nfu_rpkm_185_samples.txt.gz"), header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

cat("Non-NA/non-blank danrer_external_gene_id count:", sum(!is.na(rpkm$danrer_external_gene_id) & rpkm$danrer_external_gene_id != ""), "out of", nrow(rpkm), "\n")
cat("\nSample of actual gene symbol values:\n")
print(head(rpkm$danrer_external_gene_id[rpkm$danrer_external_gene_id != ""], 20))

# Check against our actual gene sets
nuclear_oxphos <- read_csv("data/raw/mitocarta/nuclear_oxphos_genes.csv")$gene_symbol
mtdna_genes <- c("MT-ND1","MT-ND2","MT-ND3","MT-ND4","MT-ND4L","MT-ND5","MT-ND6",
                  "MT-CO1","MT-CO2","MT-CO3","MT-CYB","MT-ATP6","MT-ATP8")
our_genes <- c(nuclear_oxphos, mtdna_genes)

# Zebrafish symbols are lowercase by convention
danio_syms_upper <- toupper(rpkm$danrer_external_gene_id)
matched <- our_genes[toupper(our_genes) %in% danio_syms_upper]
cat("\nOur genes matched in this dataset:", length(matched), "out of", length(our_genes), "\n")
cat("Matched:", paste(head(matched, 20), collapse=", "), "...\n")