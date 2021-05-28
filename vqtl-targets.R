library("data.table")
library("stringr")
set.seed(123)

# get coloc vQTLs
d <- fread("data/ldl_direct.30780.0.0.coloc.txt")

# select vQTL loci with shared casual variant with gene expression in blood
d <- d[d$PP.H4.abf > .7]
genes <- str_split(d$gene, "-", simplify=T)[,3]

# query opentargets for genes
