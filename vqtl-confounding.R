load("data/pheno.RData")

library("varGWASR")
library("dplyr")
library("data.table")
library('optparse')
source("funs.R")
set.seed(12)

option_list = list(
  make_option(c("-t", "--trait"), type="character", default=NULL, help="Variable name for outcome", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# load linker
linker <- get_filtered_linker(drop_standard_excl=TRUE, drop_non_white_british=TRUE, drop_related=TRUE, application="15825")

# load covariates
covariates <- get_covariates()
covariates$chip <- as.numeric(as.factor(covariates$chip)) - 1
pc <- get_genetic_principal_components()

# merge data
dat <- merge(linker, covariates, "appieu")
dat <- merge(dat, pheno, by.x="app15825", by.y="eid")
dat <- merge(dat, pc, "appieu")

# number of SNPs to subsample
n_snps <- 10000

# subsample genotypes
chr22 <- fread("/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/snp-stats/data.chr22.snp-stats")
chr22_s <- chr22[sample(nrow(chr22), n_snps), ]

# test for vQTL effect w/wo genetic PCs in the second-stage model
results <- data.frame()
covar <- c("age_at_recruitment.21022.0.0","sex.31.0.0", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")
for (i in 1:n_snps){
    g <- extract_variant_from_bgen(as.character(chr22_s$chromosome[i]), as.double(chr22_s$position[i]), chr22_s$alleleA[i], chr22_s$alleleB[i])
    dx <- merge(dat, g, "appieu")
    dx <- dx %>% dplyr::select(paste0("chr", chr22_s$chromosome[i], "_", chr22_s$position[i], "_", chr22_s$alleleA[i], "_", chr22_s$alleleB[i]), "age_at_recruitment.21022.0.0","sex.31.0.0", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "standing_height.50.0.0") %>% tidyr::drop_na()
    fit0 <- varGWASR::model(dx, paste0("chr", chr22_s$chromosome[i], "_", chr22_s$position[i], "_", chr22_s$alleleA[i], "_", chr22_s$alleleB[i]), "standing_height.50.0.0", covar1 = covar, covar2 = NULL)
    names(fit0) <- paste0(names(fit0), ".fit0")
    fit1 <- varGWASR::model(dx, paste0("chr", chr22_s$chromosome[i], "_", chr22_s$position[i], "_", chr22_s$alleleA[i], "_", chr22_s$alleleB[i]), "standing_height.50.0.0", covar1 = covar, covar2 = covar)
    names(fit1) <- paste0(names(fit1), ".fit1")
    results <- rbind(results, cbind(fit0, fit1))
}

# save results
write.csv(results, file="data/height-vqtl-confounding.csv")