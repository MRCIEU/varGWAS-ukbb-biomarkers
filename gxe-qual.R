load("data/pheno.RData")
library('optparse')
library('data.table')
library('dplyr')
library('broom')
library('ieugwasr')
library('robustbase')
source("funs.R")
set.seed(1234)
options(ieugwasr_api="http://64.227.44.193:8006/")

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

# split env into binary groups
dat$smoking_status.20116.0.0_b <- dat$smoking_status.20116.0.0 > 0
dat$summed_minutes_activity.22034.0.0_b <- dat$summed_minutes_activity.22034.0.0 < median(dat$summed_minutes_activity.22034.0.0, na.rm=T)
dat$alcohol_intake_frequency.1558.0.0_b <- dat$alcohol_intake_frequency.1558.0.0 < 4
dat$estimated_fat_yesterday.100004.0.0_b <- dat$estimated_fat_yesterday.100004.0.0 > median(dat$estimated_fat_yesterday.100004.0.0, na.rm=T)
dat$estimated_total_sugars_yesterday.100008.0.0_b <- dat$estimated_total_sugars_yesterday.100008.0.0 > median(dat$estimated_total_sugars_yesterday.100008.0.0, na.rm=T)
dat$age_at_recruitment.21022.0.0_b <- dat$age_at_recruitment.21022.0.0 > median(dat$age_at_recruitment.21022.0.0, na.rm=T)
dat$body_mass_index.21001.0.0_b <- dat$body_mass_index.21001.0.0 > median(dat$body_mass_index.21001.0.0, na.rm=T)
dat$sex.31.0.0_b <- dat$sex.31.0.0 > 0

# SD outcome
dat[[opt$trait]] <- dat[[opt$trait]] / sd(dat[[opt$trait]], na.rm=T)

# read in clumped vQTLs
snps <- fread(paste0("data/", opt$trait, ".clump.txt"))
stopifnot(nrow(snps)>0)

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
for (i in env_exp){
  for (j in 1:length(vqtls)){
    f <- as.formula(paste0(opt$trait, " ~ ", vqtls[j], " + age_at_recruitment.21022.0.0 + sex.31.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))
    k <- paste0(i, "_b")
    t <- tidy(lmrob(f, data=dat %>% dplyr::filter(!!sym(k) == T)))
    res1 <- t[2,]
    t <- tidy(lmrob(f, data=dat %>% dplyr::filter(!!sym(k) == F)))
    res2 <- t[2,]
    names(res1) <- paste0(names(res1), ".T")
    names(res2) <- paste0(names(res2), ".F")
    res <- cbind(res1, res2)
    results <- rbind(results, res)
  }
}

# save
write.table(results, sep="\t", quote=F, row.names=F, file=paste0("data/", opt$trait, ".gxe-qual.txt"))