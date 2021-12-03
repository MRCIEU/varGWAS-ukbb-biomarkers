library("data.table")
library("stringr")
set.seed(13)

map_snps <- function(row){
    snps <- grep("chr", stringr::str_split(row, " ")[[1]], value=T)
    return(paste0(snps, collapse="|"))
}

# gxg
gxg <- fread("data/gxg-add-finemap.txt", select=c("term", "trait", "formula"))
gxg$snps <- sapply(gxg$formula, function(row) map_snps(row))
gxg$formula <- NULL

# gxe
gxe <- fread("data/gxe-finemap.txt", select=c("term", "trait", "formula"))
gxe$snps <- sapply(gxe$formula, function(row) map_snps(row))
gxe$formula <- NULL

# all
all <- rbind(gxg, gxe)
write.csv(all, file="Table S3.csv", quote=F, row.names=F)