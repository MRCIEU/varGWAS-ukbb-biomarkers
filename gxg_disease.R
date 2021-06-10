library("data.table")
library("dplyr")
library("stringr")
source("funs.R")
set.seed(123)

# args
opt <- data.frame(trait="aspartate_aminotransferase.30650.0.0", stringsAsFactors=F)

# read in GxG associations
snps <- fread(paste0("data/", opt$trait, ".gxg.txt"))

# filter on P value
snps <- snps %>% filter(p.value < 0.05 / (1e+6 + 250000) / 30)

# split term
snps <- cbind(snps, str_split(snps$term, ":", simplify=T), stringsAsFactors=F)
v1 <- as.data.frame(str_split(snps$V1, "_", simplify=T), stringsAsFactors=F)
names(v1) <- c("chr.1", "pos.1", "oa.1", "ea.1")
v2 <- as.data.frame(str_split(snps$V2, "_", simplify=T), stringsAsFactors=F)
names(v2) <- c("chr.2", "pos.2", "oa.2", "ea.2")
snps <- cbind(snps, v1, v2, stringsAsFactors=F)
snps$chr.1 <- gsub("chr", "", snps$chr.1)
snps$chr.2 <- gsub("chr", "", snps$chr.2)

#