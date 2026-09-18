library(GEOquery)

gse <- getGEO("GSE66712", GSEMatrix = TRUE)
print(length(gse))
eset <- gse[[1]]
print(dim(eset))
print(head(pData(eset)[, c("title", "characteristics_ch1", "characteristics_ch1.1", "characteristics_ch1.2")]))

# Also list actual supplementary files available (raw counts, if any)
supp <- getGEOSuppFiles("GSE66712", makeDirectory = FALSE, fetch_files = FALSE)
print(supp)