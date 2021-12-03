load("data/pheno.RData")

library('optparse')
library('data.table')
library('dplyr')
library('broom')
library('lmtest')
library('sandwich')
library('ieugwasr')
library("varGWASR")
source("funs.R")
set.seed(1234)
options(ieugwasr_api="http://64.227.44.193:8006/")

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

# SD scale outcomes
for (e in biomarkers){
  dat[[e]] <- dat[[e]] / sd(dat[[e]], na.rm=T)
}

# read in GxE concentrating effects
d <- fread("data/gxe-qual-conc.txt", select=c("key"))
d <- cbind(d, as.data.frame(str_split(d$key, ":", simplify=T), stringsAsFactors=F) %>% dplyr::rename(u="V1", x="V2", y="V3"))

# drop conc effects that attentuate when adjusted for main effects (see check_gxe_sub.R)
d <- d %>% dplyr::filter(key == "body_mass_index.21001.0.0:chr9_132566666_G_A:alanine_aminotransferase.30620.0.0" | key == "body_mass_index.21001.0.0:chr22_44324727_C_G:triglycerides.30870.0.0")
d$key <- NULL

# add in two large quantiative effects
d <- rbind(d, data.frame(
  u="body_mass_index.21001.0.0",
  x="chr22_44324727_C_G",
  y="alanine_aminotransferase.30620.0.0",
  stringsAsFactors=F
))
d <- rbind(d, data.frame(
  u="sex.31.0.0",
  x="chr4_9926051_A_G",
  y="urate.30880.0.0",
  stringsAsFactors=F
))
d$key <- paste0(d$u, ":", d$x, ":", d$y)

# read in fine-mapped covarites
fm <- fread("Table S3.csv")
fm$key <- paste0(fm$term, ":", fm$trait)
names(fm)[3] <- "covar"
fm <- fm %>% dplyr::select(key, covar)
d <- merge(d, fm, "key")

# create snp list
snps <- d$x
snps <- c(snps, stringr::str_split(d$covar, "\\|") %>% unlist)
snps <- as.data.frame(str_split(snps, "_", simplify=T), stringsAsFactors=F) %>% dplyr::rename(chr="V1", pos="V2", oa="V3", ea="V4")
snps$chr <- gsub("chr", "", snps$chr)
snps <- unique(snps)

# load dosages
for (i in 1:nrow(snps)){
    dosage <- extract_variant_from_bgen(as.character(snps$chr[i]), as.double(snps$pos[i]), snps$oa[i], snps$ea[i])
    dat <- merge(dat, dosage, "appieu")
}

# test for effect of SNP stratified by modifier
results <- data.frame()
results_var <- data.frame()
for (i in 1:nrow(d)){
  # define modifier group
  k <- paste0(d$u[i], "_b")

  # define linear model
  covar_snps <- str_split(d$covar[i], "\\|")[[1]]
  covar_snps <- covar_snps[d$x[i] != covar_snps]
  if (length(covar_snps) > 0){
    f <- as.formula(paste0(d$y[i], " ~ ", d$x[i], " + age_at_recruitment.21022.0.0 + sex.31.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + ", paste0(covar_snps, collapse="+")))
  } else {
    f <- as.formula(paste0(d$y[i], " ~ ", d$x[i], " + age_at_recruitment.21022.0.0 + sex.31.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))
  }

  # split dataset
  dat0 <- dat %>% dplyr::filter(!!sym(k) == T)
  dat1 <- dat %>% dplyr::filter(!!sym(k) == F)

  # test SNP effect on outcome by stratified modifier
  mod <- lm(f, data=dat0)
  t <- coeftest(mod, vcov = vcovHC(mod, type = "HC0")) %>% tidy
  t <- t[2,]
  t$mod <- 0
  t$mod_pheno <- k
  t$trait <- d$y[i]
  results <- rbind(results, t)

  mod <- lm(f, data=dat1)
  t <- coeftest(mod, vcov = vcovHC(mod, type = "HC0")) %>% tidy
  t <- t[2,]
  t$mod <- 1
  t$mod_pheno <- k
  t$trait <- d$y[i]
  results <- rbind(results, t)

  # test SNP effect on outcome variance ajusted for int
  dat$XU <- dat[[d$x[i]]] * dat[[d$u[i]]]
  covar <- c(
    "age_at_recruitment.21022.0.0",
    "sex.31.0.0",
    "PC1",
    "PC2",
    "PC3",
    "PC4",
    "PC5",
    "PC6",
    "PC7",
    "PC8",
    "PC9",
    "PC10",
    covar_snps
  )
  covar <- unique(covar)
  dat2 <- dat %>% dplyr::select(!!d$x[i], !!covar, !!d$y[i]) %>% tidyr::drop_na()
  fit <- model(dat2, d$x[i], d$y[i], covar1 = covar, covar2 = covar)
  fit$int=F
  fit$term <- paste0(d$x[i], ":", d$u[i])
  fit$trait <- d$y[i]
  results_var <- rbind(results_var, fit)
  covar <- c(covar,d$u[i],"XU")
  covar <- unique(covar)
  dat2 <- dat %>% dplyr::select(!!d$x[i], !!covar, !!d$y[i]) %>% tidyr::drop_na()
  fit <- model(dat2, d$x[i], d$y[i], covar1 = covar, covar2 = covar)
  fit$trait <- d$y[i]
  fit$int=T
  fit$term <- paste0(d$x[i], ":", d$u[i])
  results_var <- rbind(results_var, fit)
}

# save
write.table(results, sep="\t", quote=F, row.names=F, file=paste0("data/gxe-qual-var1.txt"))
write.table(results_var, sep="\t", quote=F, row.names=F, file=paste0("data/gxe-qual-var2.txt"))