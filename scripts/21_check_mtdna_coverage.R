rpkm <- read.table(gzfile("data/raw/GSE66712_nfu_rpkm_185_samples.txt.gz"), header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

mtdna_pattern <- "^mt-|^nd[1-6]$|^nd4l$|^co[1-3]$|^cytb$|^atp6$|^atp8$"
hits <- rpkm$danrer_external_gene_id[grepl(mtdna_pattern, rpkm$danrer_external_gene_id, ignore.case = TRUE)]
cat("Possible mtDNA gene matches:\n")
print(unique(hits))