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
snps$key1 <- paste0(snps$chr.1, ":", snps$pos.1)
snps$key2 <- paste0(snps$chr.2, ":", snps$pos.2)

# load list of genes for these variants
genes <- fread("data/aspartate_aminotransferase.30650.0.0.gxg-usnps-closest.txt", select=c("V1", "V3", "V7", "V8"), col.names=c("chr", "pos", "ensg", "gene"))
genes$key <- paste0(genes$chr, ":", genes$pos)
genes$chr <- NULL
genes$pos <- NULL
genes$ensg <- NULL

# merge
snps <- merge(snps, genes, by.x="key1", by.y="key")
snps <- merge(snps, genes, by.x="key2", by.y="key", suffixes=c(".1", ".2"))

# add gxg
snps$gxg.x <- paste0(snps$gene.1, "|", snps$gene.2)
snps$gxg.y <- paste0(snps$gene.2, "|", snps$gene.1)

snps %>% select(gxg)