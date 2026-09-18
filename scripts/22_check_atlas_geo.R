library(GEOquery)

supp <- getGEOSuppFiles("GSE308970", makeDirectory = FALSE, fetch_files = FALSE)
print(supp)