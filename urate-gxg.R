library("ieugwasr")
library('optparse')
library('data.table')
library('dplyr')
library('broom')
source("funs.R")
set.seed(1234)

# get tophits for urate
urate <- tophits("ukb-d-30880_irnt")

# subset two vQTL regions
pkd2_abcg2 <- urate %>% filter(chr==4 & position >= 89050026-500000 & position <= 89050026+500000)
slc2a9_wdr1 <- urate %>% filter(chr==4 & position >= 9926051-500000 & position <= 9926051+500000)

# read in extracted phenotypes
pheno <- fread("data/urate.30880.0.0.txt")

# load dosages
snps <- rbind(pkd2_abcg2, slc2a9_wdr1)
for (i in 1:nrow(snps)){
    dosage <- extract_variant_from_bgen(as.character(snps$chr[i]), as.double(snps$position[i]), snps$nea[i], snps$ea[i])
    pheno <- merge(pheno, dosage, "appieu")
}

# select vQTLs
variants <- grep("^chr", names(pheno), value=T)

# test for interaction between each snp
results <- data.frame()
for (i in 1:length(variants)){
  for (j in 1:length(variants)){
    f <- as.formula(paste0("urate.30880.0.0 ~ age_at_recruitment.21022.0.0 + sex.31.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + ", paste0(variants[i], " * " ,variants[j], collapse=" + ")))
    fit <- lm(f, pheno)
    t <- tidy(fit)
    results <- rbind(results, t[grep(":", t$term),])
  }
}