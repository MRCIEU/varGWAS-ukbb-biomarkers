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
pc <- get_genetic_principal_components()

# merge data
dat <- merge(linker, covariates, "appieu")
dat <- merge(dat, pheno, by.x="app15825", by.y="eid")
dat <- merge(dat, pc, "appieu")

# select fields for GWAS
dat <- dat[,c("appieu", "sex.31.0.0", "age_at_recruitment.21022.0.0", opt$trait, paste0("PC", seq(1, 10))), with=F]

# drop missing values
dat <- dat[complete.cases(dat), ]

# SD scale
dat[[opt$trait]] <- dat[[opt$trait]] / sd(dat[[opt$trait]])
dat$age_at_recruitment.21022.0.0 <- dat$age_at_recruitment.21022.0.0 / sd(dat$age_at_recruitment.21022.0.0)

# write out pheno for vGWAS
write.table(dat, file=paste0("data/", opt$trait, ".txt", row.names=F, quote=F, sep=",")