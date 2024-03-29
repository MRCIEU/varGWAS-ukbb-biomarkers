load("data/pheno.RData")

library('optparse')
library('data.table')
library('dplyr')
library('car')
library('broom')
library('ieugwasr')
library("lmtest")
library("sandwich")
library("varGWASR")
source("funs.R")
set.seed(1234)
options(ieugwasr_api="http://web-dc1-bms-d0.infra.bris.ac.uk:5002/")

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

# read in GxG effects
d <- fread("data/gxg.txt")
d <- d %>% dplyr::filter(trait != "body_mass_index.21001.0.0")
d <- d %>% dplyr::filter(trait != "c_reactive_protein.30710.0.0") # not replicating on log scale
d <- d %>% dplyr::filter(trait != "triglycerides.30870.0.0") # not replicating on log scale
d <- d %>% dplyr::filter(p.value < 5e-8)
d <- cbind(d, as.data.frame(str_split(d$term, ":", simplify=T), stringsAsFactors=F))
d$V3 <- d$V1
d$V1 <- d$V2
d$V2 <- d$V3
d$V3 <- NULL
v1 <- as.data.frame(str_split(d$V1, "_", simplify=T), stringsAsFactors=F)
names(v1) <- c("chr.1", "pos.1", "oa.1", "ea.1")
v2 <- as.data.frame(str_split(d$V2, "_", simplify=T), stringsAsFactors=F)
names(v2) <- c("chr.2", "pos.2", "oa.2", "ea.2")
d <- cbind(d, v1, v2)
d$chr.1 <- gsub("chr", "", d$chr.1)
d$chr.2 <- gsub("chr", "", d$chr.2)

# load dosages
snps <- d %>% dplyr::select(chr.1, pos.1, oa.1, ea.1) %>% rename(chr="chr.1", pos="pos.1", oa="oa.1", ea="ea.1")
snps <- rbind(snps, d %>% dplyr::select(chr.2, pos.2, oa.2, ea.2) %>% rename(chr="chr.2", pos="pos.2", oa="oa.2", ea="ea.2"))
snps <- unique(snps)
snps$key <- paste0("chr", snps$chr, "_", snps$pos, "_", snps$oa, "_", snps$ea)

for (i in 1:nrow(snps)){
    dosage <- extract_variant_from_bgen(as.character(snps$chr[i]), as.double(snps$pos[i]), snps$oa[i], snps$ea[i])
    dat <- merge(dat, dosage, "appieu")
}

# test for effect of SNP stratified by modifier
results <- data.frame()
results_var <- data.frame()
for (i in 1:nrow(d)){
  # set Z
  s <- sd(dat[[d$trait[i]]], na.rm=T)
  dat$z <- (dat[[d$trait[i]]] - mean (dat[[d$trait[i]]], na.rm=T)) / s
  dat$z <- abs(dat$z)

  # SD scale
  dat$out <- dat[[d$trait[i]]] / s

  # modifier
  k <- d$V2[i]
  dat$mod_gt <- round(dat[[k]])
  dat0 <- dat %>% dplyr::filter(mod_gt == 0)
  dat1 <- dat %>% dplyr::filter(mod_gt == 1)
  dat2 <- dat %>% dplyr::filter(mod_gt == 2)

  # set fit
  f <- as.formula(paste0("out ~ ", d$V1[i], " + age_at_recruitment.21022.0.0 + sex.31.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))

  # test SNP effect on outcome by stratified modifier
  mod <- lm(f, data=dat0)
  t <- coeftest(mod, vcov = vcovHC(mod, type = "HC0")) %>% tidy
  t <- t[2,]
  t$mod <- 0
  t$mod_snp <- k
  t$trait <- d$trait[i]
  results <- rbind(results, t)

  mod <- lm(f, data=dat1)
  t <- coeftest(mod, vcov = vcovHC(mod, type = "HC0")) %>% tidy
  t <- t[2,]
  t$mod <- 1
  t$mod_snp <- k
  t$trait <- d$trait[i]
  results <- rbind(results, t)

  mod <- lm(f, data=dat2)
  t <- coeftest(mod, vcov = vcovHC(mod, type = "HC0")) %>% tidy
  t <- t[2,]
  t$mod <- 2
  t$mod_snp <- k
  t$trait <- d$trait[i]
  results <- rbind(results, t)

  # test SNP effect on outcome variance w/wo adjusting for interaction
  dat$XU <- dat[[d$V1[i]]] * dat[[d$V2[i]]]
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
    "PC10"
  )
  # without interaction
  dat2 <- dat %>% dplyr::filter(z < 5) %>% dplyr::select(!!d$V1[i], !!covar, out) %>% tidyr::drop_na()
  fit <- varGWASR::model(dat2, d$V1[i], "out", covar1 = covar, covar2 = covar)
  fit$int=F
  fit$term <- paste0(d$V1[i], ":", d$V2[i])
  fit$trait <- d$trait[i]
  results_var <- rbind(results_var, fit)
  
  # with interaction
  covar <- c(covar,d$V2[i],"XU")
  dat2 <- dat %>% dplyr::filter(z < 5) %>% dplyr::select(!!d$V1[i], !!covar, out) %>% tidyr::drop_na()
  fit <- varGWASR::model(dat2, d$V1[i], "out", covar1 = covar, covar2 = covar)
  fit$trait <- d$trait[i]
  fit$int=T
  fit$term <- paste0(d$V1[i], ":", d$V2[i])
  results_var <- rbind(results_var, fit)
}

# save
write.table(results, sep="\t", quote=F, row.names=F, file=paste0("data/gxg-qual.txt"))
write.table(results_var, sep="\t", quote=F, row.names=F, file=paste0("data/gxg-qual-var.txt"))