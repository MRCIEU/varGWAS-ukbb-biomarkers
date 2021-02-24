load("data/pheno.RData")
library('data.table')
library('dplyr')
source("funs.R")
set.seed(1234)

# export phenotypes for vGWAS

# load linker
linker <- get_filtered_linker(drop_standard_excl=TRUE, drop_non_white_british=TRUE, drop_related=TRUE, application="15825")

# load covariates
covariates <- get_covariates()

# merge data
dat <- merge(linker, covariates, "appieu")
dat <- merge(dat, pheno, by.x="app15825", by.y="eid")

# select fields for GWAS
dat <- dat[,c("appieu", "appieu", "sex.31.0.0", "chip", "age_at_recruitment.21022.0.0", "body_mass_index.21001.0.0"), with=F]
names(dat)[1] <- "FID"
names(dat)[2] <- "IID"

# drop missing values
dat <- dat[complete.cases(dat), ]

# SD scale
dat$body_mass_index.21001.0.0 <- dat$body_mass_index.21001.0.0 / sd(dat$body_mass_index.21001.0.0)

# write out pheno for vGWAS
write.table(dat, file="data/bmi.txt", row.names=F, quote=F, sep=",")