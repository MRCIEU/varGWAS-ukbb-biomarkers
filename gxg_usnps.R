library("data.table")
library("dplyr")
library("stringr")
source("funs.R")
set.seed(123)

# args
opt <- data.frame(trait="aspartate_aminotransferase.30650.0.0", stringsAsFactors=F)

# load all vGWAS for trait
vgwas <- get_variants(opt$trait)
vgwas$key <- paste0("chr", vgwas$key)

# read in GxG associations
snps <- fread(paste0("data/", opt$trait, ".gxg.txt"))

# filter on P value
snps <- snps %>% filter(p.value < 0.05 / (1e+6 + 250000) / 30)

# split term
snps <- cbind(snps, str_split(snps$term, ":", simplify=T), stringsAsFactors=F)
usnps <- unique(c(snps$V1,snps$V2), stringsAsFactors=F)

# subset vGWAS for these variants
vgwas <- vgwas %>% filter(key %in% usnps)

# write to file for coloc
write.table(vgwas, file=paste0("data/", opt$trait, ".gxg-usnps.txt"), sep="\t", quote=F, row.name=F)