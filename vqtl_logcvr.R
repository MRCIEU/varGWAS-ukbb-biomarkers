library("ieugwasr")
library("dplyr")
library("multcomp")
library("broom")
source("funs.R")
set.seed(1234)
options(ieugwasr_api="http://64.227.44.193:8006/")

bp <- function(dat, snp, outcome){
    dat$x <- dat[[snp]]
    dat$xsq <- dat$x^2
    fit1 <- lm(paste0(outcome, " ~ x + age_at_recruitment.21022.0.0 + sex.31.0.0 + ", paste0("PC", seq(1,10), collapse= " + ")), data=dat)
    dat$dsq <- resid(fit1)^2
    fit2 <- lm(dsq ~ x + xsq, data=dat)
    beta0 <- glht(model=fit1, linfct=paste("Intercept == 0"))
    beta1 <- glht(model=fit1, linfct=paste("Intercept + x*1 == 0"))
    beta2 <- glht(model=fit1, linfct=paste("Intercept + x*2 == 0"))
    varbeta0 <- glht(model=fit2, linfct=paste("Intercept == 0"))
    varbeta1 <- glht(model=fit2, linfct=paste("Intercept + x*1 + xsq*1 == 0"))
    varbeta2 <- glht(model=fit2, linfct=paste("Intercept + x*2 + xsq*4 == 0"))
    res <- cbind(
        tidy(beta0) %>% dplyr::select("estimate", "std.error") %>% rename(estimate.0="estimate", std.error.0="std.error"),
        tidy(beta1) %>% dplyr::select("estimate", "std.error") %>% rename(estimate.1="estimate", std.error.1="std.error"),
        tidy(beta2) %>% dplyr::select("estimate", "std.error") %>% rename(estimate.2="estimate", std.error.2="std.error"),
        tidy(varbeta0) %>% dplyr::select("estimate", "std.error") %>% rename(var.estimate.0="estimate", var.std.error.0="std.error"),
        tidy(varbeta1) %>% dplyr::select("estimate", "std.error") %>% rename(var.estimate.1="estimate", var.std.error.1="std.error"),
        tidy(varbeta2) %>% dplyr::select("estimate", "std.error") %>% rename(var.estimate.2="estimate", var.std.error.2="std.error"),
        n0=table(round(dat$x))[1],
        n1=table(round(dat$x))[2],
        n2=table(round(dat$x))[3]
    )
    res$sd.estimate.0 <- sqrt(res$var.estimate.0)
    res$sd.estimate.1 <- sqrt(res$var.estimate.1)
    res$sd.estimate.2 <- sqrt(res$var.estimate.2)
    res$snp <- snp
    res$outcome <- outcome
    return(res)
}

# load phenotypes
load("data/pheno.RData")
linker <- get_filtered_linker(drop_standard_excl=TRUE, drop_non_white_british=TRUE, drop_related=TRUE, application="15825")
covariates <- get_covariates()
pc <- get_genetic_principal_components()
dat <- merge(linker, covariates, "appieu")
dat <- merge(dat, pheno, by.x="app15825", by.y="eid")
dat <- merge(dat, pc, "appieu")
ldl.pheno <- dat %>% dplyr::select(appieu, age_at_recruitment.21022.0.0, sex.31.0.0, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, ldl_direct.30780.0.0)
ldl.pheno <- na.omit(ldl.pheno)
hdl.pheno <- dat %>% dplyr::select(appieu, age_at_recruitment.21022.0.0, sex.31.0.0, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, hdl_cholesterol.30760.0.0)
hdl.pheno <- na.omit(hdl.pheno)

# load instruments for lipids at drug target loci & cross ref with vGWAS

# LDL-c

# get instruments from non-UKBB
ldl <- tophits("ieu-a-300", clump=0)

# dplyr::filter drug target loci
hmgcr <- ldl %>% dplyr::filter(chr == "5" & position > 74632154-500000 & position < 74657929+500000 & eaf >= 0.1 & eaf <= 0.9) %>% dplyr::select(rsid, p) %>% rename(pval=p) %>% ld_clump
hmgcr <- ldl %>% dplyr::filter(rsid %in% hmgcr$rsid) %>% mutate(gene="HMGCR")
pcsk9 <- ldl %>% dplyr::filter(chr == "1" & position > 55505221-500000 & position < 55530525+500000 & eaf >= 0.1 & eaf <= 0.9) %>% dplyr::select(rsid, p) %>% rename(pval=p) %>% ld_clump
pcsk9 <- ldl %>% dplyr::filter(rsid %in% pcsk9$rsid) %>% mutate(gene="PCSK9")
npc1l1 <- ldl %>% dplyr::filter(chr == "7" & position > 44552134-500000 & position < 44580914+500000 & eaf >= 0.1 & eaf <= 0.9) %>% dplyr::select(rsid, p) %>% rename(pval=p) %>% ld_clump
npc1l1 <- ldl %>% dplyr::filter(rsid %in% npc1l1$rsid) %>% mutate(gene="NPC1L1")
cetp <- ldl %>% dplyr::filter(chr == "16" & position > 56995762-500000 & position < 57017757+500000 & eaf >= 0.1 & eaf <= 0.9) %>% dplyr::select(rsid, p) %>% rename(pval=p) %>% ld_clump
cetp <- ldl %>% dplyr::filter(rsid %in% cetp$rsid) %>% mutate(gene="CETP")

# orient effects to be lipid lowering
ldl <- rbind(hmgcr, pcsk9, npc1l1, cetp)
ldl$flip <- ldl$beta > 0
ldl$allele <- ""
ldl[which(ldl$flip),]$beta <- ldl[which(ldl$flip),]$beta * -1
ldl[which(ldl$flip),]$eaf <- 1 - ldl[which(ldl$flip),]$eaf
ldl[which(ldl$flip),]$allele <- ldl[which(ldl$flip),]$nea
ldl[which(ldl$flip),]$nea <- ldl[which(ldl$flip),]$ea
ldl[which(ldl$flip),]$ea <- ldl[which(ldl$flip),]$allele

# test for SNP effect on LDL-c variance
for (i in 1:nrow(ldl)){
    dosage <- extract_variant_from_bgen(as.character(ldl$chr[i]), as.double(ldl$position[i]), ldl$nea[i], ldl$ea[i])
    names(dosage)[1] <- ldl$rsid[i]
    ldl.pheno <- merge(ldl.pheno, dosage, "appieu")
}

ldl.results <- data.frame()
for (i in 1:nrow(ldl)){
    ldl.results <- rbind(ldl.results, bp(ldl.pheno, ldl$rsid[i], "ldl_direct.30780.0.0"))
}

# HDL-c

# get instruments
hdl <- tophits("ieu-a-299", clump=0)

# dplyr::filter drug target loci
cetp <- hdl %>% dplyr::filter(chr == "16" & position > 56995762-500000 & position < 57017757+500000 & eaf >= 0.1 & eaf <= 0.9) %>% dplyr::select(rsid, p) %>% rename(pval=p) %>% ld_clump
cetp <- hdl %>% dplyr::filter(rsid %in% cetp$rsid) %>% mutate(gene="CETP")

# orient effects to be lipid lowering
hdl <- cetp
hdl$flip <- hdl$beta < 0
hdl$allele <- ""
hdl[which(hdl$flip),]$beta <- hdl[which(hdl$flip),]$beta * -1
hdl[which(hdl$flip),]$eaf <- 1 - hdl[which(hdl$flip),]$eaf
hdl[which(hdl$flip),]$allele <- hdl[which(hdl$flip),]$nea
hdl[which(hdl$flip),]$nea <- hdl[which(hdl$flip),]$ea
hdl[which(hdl$flip),]$ea <- hdl[which(hdl$flip),]$allele

# test for SNP effect on HDL-c variance
pheno <- fread("data/hdl_cholesterol.30760.0.0.txt")
for (i in 1:nrow(hdl)){
    dosage <- extract_variant_from_bgen(as.character(hdl$chr[i]), as.double(hdl$position[i]), hdl$nea[i], hdl$ea[i])
    names(dosage)[1] <- hdl$rsid[i]
    pheno <- merge(pheno, dosage, "appieu")
}

hdl.results <- data.frame()
for (i in 1:nrow(hdl)){
    hdl.results <- rbind(hdl.results, bp(pheno, hdl$rsid[i], "hdl_cholesterol.30760.0.0"))
}

# write results
write.table(file="vqtl_logcvr.txt", rbind(ldl.results, hdl.results), sep="\t", row.names=F, quote=F)