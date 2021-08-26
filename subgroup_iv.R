library("ieugwasr")
library("dplyr")
library("multcomp")
library("broom")
library('robustbase')
source("funs.R")
set.seed(1234)
options(ieugwasr_api="http://64.227.44.193:8006/")

# load phenotypes
load("data/pheno.RData")
linker <- get_filtered_linker(drop_standard_excl=TRUE, drop_non_white_british=TRUE, drop_related=TRUE, application="15825")
covariates <- get_covariates()
pc <- get_genetic_principal_components()
dat <- merge(linker, covariates, "appieu")
dat <- merge(dat, pheno, by.x="app15825", by.y="eid")
dat <- merge(dat, pc, "appieu")

# get eQTL instruments
lpl <- fread("data/lpl.txt")
lpl <- lpl %>% dplyr::filter(chr == "8" & position > 19759228-500000 & position < 19824769+500000)
lpl$gene <- "LPL"

slc2a9 <- fread("data/slc2a9.txt")
slc2a9 <- slc2a9 %>% dplyr::filter(chr == "4" & position > 9772777-500000 & position < 10056560+500000)
slc2a9$gene <- "SLC2A9"

# combine
snps <- rbind(lpl, slc2a9)

results <- data.frame()
for (i in 1:nrow(snps)){
    # load dosage
    dosage <- extract_variant_from_bgen(as.character(snps$chr[i]), as.double(snps$position[i]), snps$nea[i], snps$ea[i])
    names(dosage)[1] <- snps$rsid[i]
    dat <- merge(dat, dosage, "appieu")

    # test for SNP-outcome effect stratified by modifier
    f <- as.formula(paste0(snps$rsid[i], " + age_at_recruitment.21022.0.0 + sex.31.0.0 + PC1 + PC2 + PC3 + PC4 + PC6 + PC7 + PC8 + PC9 + PC10"))
    fit <- tidy(lmrob(f, data=dat))[2]
}