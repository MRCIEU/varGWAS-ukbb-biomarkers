load("data/pheno.RData")
library('data.table')
library('dplyr')
library('optparse')
source("funs.R")
set.seed(1234)

option_list = list(
  make_option(c("-t", "--trait"), type="character", default=NULL, help="Variable name for outcome", metavar="character")
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

# select fields for GWAS
dat <- dat[,c("appieu", "sex.31.0.0", "age_at_recruitment.21022.0.0", opt$trait, paste0("PC", seq(1, 10))), with=F]

# drop missing values
dat <- dat[complete.cases(dat), ]

# drop outlier values
s <- sd(dat[[opt$trait]])
dat$z <- (dat[[opt$trait]] - mean (dat[[opt$trait]])) / s
dat$z <- abs(dat$z)
dat <- dat %>% dplyr::filter(z < 5)
dat$z <- NULL

# SD scale
dat[[opt$trait]] <- dat[[opt$trait]] / s

# squared terms
dat$age_at_recruitment.21022.0.0_sq <- dat$age_at_recruitment.21022.0.0^2
dat$PC1_sq <- dat$PC1^2
dat$PC2_sq <- dat$PC2^2
dat$PC3_sq <- dat$PC3^2
dat$PC4_sq <- dat$PC4^2
dat$PC5_sq <- dat$PC5^2
dat$PC6_sq <- dat$PC6^2
dat$PC7_sq <- dat$PC7^2
dat$PC8_sq <- dat$PC8^2
dat$PC9_sq <- dat$PC9^2
dat$PC10_sq <- dat$PC10^2

# write out pheno for vGWAS
write.table(dat, file=paste0("data/", opt$trait, ".txt"), row.names=F, quote=F, sep=",")