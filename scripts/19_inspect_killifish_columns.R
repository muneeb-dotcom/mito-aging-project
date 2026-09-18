rpkm <- read.table(gzfile("data/raw/GSE66712_nfu_rpkm_185_samples.txt.gz"), header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
cat("Full column list:\n")
print(colnames(rpkm))
cat("\nSample orthology_mapping values:\n")
print(table(rpkm$orthology_mapping))
cat("\nExample rows for our genes of interest (search danrer_external_gene_id for nampt/sirt1):\n")
print(rpkm[grepl("nampt|sirt1|ppargc1a", rpkm$danrer_external_gene_id, ignore.case = TRUE),
           c("danrer_external_gene_id", "orthology_mapping")])