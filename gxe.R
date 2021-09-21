load("data/pheno.RData")
library('optparse')
library('data.table')
library('dplyr')
library('broom')
library('ieugwasr')
library("lmtest")
library("sandwich")
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

# SD scale env
for (e in env_exp){
  dat[[e]] <- dat[[e]] / sd(dat[[e]], na.rm=T)
}
dat[[paste0(opt$trait, "_log")]] <- log(dat[[opt$trait]])
dat[[opt$trait]] <- dat[[opt$trait]] / sd(dat[[opt$trait]], na.rm=T)
dat[[paste0(opt$trait, "_log")]] <- dat[[paste0(opt$trait, "_log")]] / sd(dat[[paste0(opt$trait, "_log")]], na.rm=T)

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
results_log <- data.frame()
for (i in env_exp){
  for (j in 1:length(vqtls)){
    # test GxE
    message("Testing GxE for: ", i, " ", vqtls[j])

    f <- as.formula(paste0(opt$trait, " ~ age_at_recruitment.21022.0.0 + sex.31.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + ", paste0(i, " * " ,vqtls[j], collapse=" + ")))
    mod <- lm(f, data=dat)
    t <- coeftest(mod, vcov = vcovHC(mod, type = "HC0")) %>% tidy
    t <- t %>% dplyr::filter(grepl(":", t$term))
    results <- rbind(results, t)

    f <- as.formula(paste0(opt$trait, "_log ~ age_at_recruitment.21022.0.0 + sex.31.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + ", paste0(i, " * " ,vqtls[j], collapse=" + ")))
    mod <- lm(f, data=dat)
    t <- coeftest(mod, vcov = vcovHC(mod, type = "HC0")) %>% tidy
    t <- t %>% dplyr::filter(grepl(":", t$term))
    results_log <- rbind(results_log, t)
  }
}

# save
write.table(results, sep="\t", quote=F, row.names=F, file=paste0("data/", opt$trait, ".gxe.txt"))
write.table(results_log, sep="\t", quote=F, row.names=F, file=paste0("data/", opt$trait, ".gxe-log.txt"))