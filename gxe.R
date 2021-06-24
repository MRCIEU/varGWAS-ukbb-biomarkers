load("data/pheno.RData")
library('optparse')
library('data.table')
library('dplyr')
library('broom')
library('ieugwasr')
source("funs.R")
set.seed(1234)

option_list = list(
  make_option(c("-t", "--trait"), type="character", default=NULL, help="Name of trait", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# export phenotypes for vGWAS
message(opt$trait)

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

# read in vGWAS
snps <- get_variants(opt$trait)

# filter on P value
vqtls <- snps %>% filter(phi_p < (5e-8/30)) %>% select("rsid", "phi_p") %>% rename(pval = phi_p)

# clump records
vqtls <- ld_clump(vqtls)
snps <- snps[snps$rsid %in% vqtls$rsid]

# add key
snps$key <- paste0("chr", snps$chr, "_", snps$pos, "_", snps$oa, "_", snps$ea)

# load dosages
for (i in 1:nrow(snps)){
    dosage <- extract_variant_from_bgen(as.character(snps$chr[i]), as.double(snps$pos[i]), snps$oa[i], snps$ea[i])
    dat <- merge(dat, dosage, "appieu")
}

# select vQTLs
vqtls <- grep("^chr", names(dat), value=T)

# test for interaction between each snp
results <- data.frame()
for (i in c("smoking_status.20116.0.0","summed_minutes_activity.22034.0.0","alcohol_intake_frequency.1558.0.0","estimated_fat_yesterday.100004.0.0","estimated_total_sugars_yesterday.100008.0.0","sex.31.0.0","age_at_recruitment.21022.0.0","body_mass_index.21001.0.0")){
  for (j in 1:length(vqtls)){
    # test GxE
    message("Testing GxE for: ", i, " ", vqtls[j])
    f <- as.formula(paste0(opt$trait, " ~ age_at_recruitment.21022.0.0 + sex.31.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + ", paste0(i, " * " ,vqtls[j], collapse=" + ")))
    fit <- lm(f, dat)
    t <- tidy(fit)

    # store results
    results <- rbind(results, t[grep(":", t$term),])
  }
}

# save
write.table(results, sep="\t", quote=F, row.names=F, file=paste0("data/", opt$trait, ".gxe.txt"))