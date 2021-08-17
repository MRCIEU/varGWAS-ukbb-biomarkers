library("ieugwasr")
library("dplyr")
library("multcomp")
library("broom")
library("tidyverse")
library("data.table")
library("IRanges")
library("stringr")
library("GenomicRanges")
library("TwoSampleMR")
library("ieugwasr")
source("funs.R")
set.seed(1234)
options(ieugwasr_api="http://64.227.44.193:8006/")

bp <- function(dat, snp, outcome){
    dat <- dat %>% dplyr::select(!!snp, !!outcome, age_at_recruitment.21022.0.0, sex.31.0.0, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10) %>% tidyr::drop_na()
    dat$x <- dat[[snp]]
    dat$xsq <- dat$x^2
    fit1 <- lm(paste0(outcome, " ~ x + age_at_recruitment.21022.0.0 + sex.31.0.0 + ", paste0("PC", seq(1,10), collapse= " + ")), data=dat)
    dat$dsq <- resid(fit1)^2
    fit2 <- lm(paste0("dsq ~ x + xsq + age_at_recruitment.21022.0.0 + sex.31.0.0 + ", paste0("PC", seq(1,10), collapse= " + ")), data=dat)
    fitnull <- lm(paste0("dsq ~ 1 + age_at_recruitment.21022.0.0 + sex.31.0.0 + ", paste0("PC", seq(1,10), collapse= " + ")), data=dat)
    beta0 <- glht(model=fit1, linfct=paste("Intercept == 0"))
    beta1 <- glht(model=fit1, linfct=paste("Intercept + x*1 == 0"))
    beta2 <- glht(model=fit1, linfct=paste("Intercept + x*2 == 0"))
    varbeta0 <- glht(model=fit2, linfct=paste("Intercept == 0"))
    varbeta1 <- glht(model=fit2, linfct=paste("Intercept + x*1 + xsq*1 == 0"))
    varbeta2 <- glht(model=fit2, linfct=paste("Intercept + x*2 + xsq*4 == 0"))
    varbeta1_est <- glht(model=fit2, linfct=paste("x*1 + xsq*1 == 0"))
    varbeta2_est <- glht(model=fit2, linfct=paste("x*2 + xsq*4 == 0"))
    res <- cbind(
        tidy(beta0) %>% dplyr::select("estimate", "std.error") %>% dplyr::rename(estimate.0="estimate", std.error.0="std.error"),
        tidy(beta1) %>% dplyr::select("estimate", "std.error") %>% dplyr::rename(estimate.1="estimate", std.error.1="std.error"),
        tidy(beta2) %>% dplyr::select("estimate", "std.error") %>% dplyr::rename(estimate.2="estimate", std.error.2="std.error"),
        tidy(varbeta0) %>% dplyr::select("estimate", "std.error") %>% dplyr::rename(var.estimate.0="estimate", var.std.error.0="std.error"),
        tidy(varbeta1) %>% dplyr::select("estimate", "std.error") %>% dplyr::rename(var.estimate.1="estimate", var.std.error.1="std.error"),
        tidy(varbeta2) %>% dplyr::select("estimate", "std.error") %>% dplyr::rename(var.estimate.2="estimate", var.std.error.2="std.error"),
        tidy(varbeta1_est) %>% dplyr::select("estimate", "std.error") %>% dplyr::rename(var_effect.estimate.1="estimate", var_effect.std.error.1="std.error"),
        tidy(varbeta2_est) %>% dplyr::select("estimate", "std.error") %>% dplyr::rename(var_effect.estimate.2="estimate", var_effect.std.error.2="std.error"),
        n0=table(round(dat$x))[1],
        n1=table(round(dat$x))[2],
        n2=table(round(dat$x))[3],
        p=tidy(anova(fitnull, fit2))$p.value[2]
    )
    res$snp <- snp
    res$outcome <- outcome
    return(res)
}

get_lci <- function(b, se){
    b <- as.numeric(b)
    se <- as.numeric(se)
    lci <- b - (1.96 * se)
    return(lci)
}
get_uci <- function(b, se){
    b <- as.numeric(b)
    se <- as.numeric(se)
    uci <- b + (1.96 * se)
    return(uci)
}

# load phenotypes
load("data/pheno.RData")
linker <- get_filtered_linker(drop_standard_excl=TRUE, drop_non_white_british=TRUE, drop_related=TRUE, application="15825")
covariates <- get_covariates()
pc <- get_genetic_principal_components()
dat <- merge(linker, covariates, "appieu")
dat <- merge(dat, pheno, by.x="app15825", by.y="eid")
dat <- merge(dat, pc, "appieu")

# extract top hits for GLGC lipid traits from OpenGWAS
d <- rbind(
    tophits("ieu-a-299"),
    tophits("ieu-a-300"),
    tophits("ieu-a-302")
)
d$key <- paste0("chr", d$chr, "_", d$position, "_", d$nea, "_", d$ea)

# subset snps in drug target loci
loci <- data.frame()
res <- d %>% dplyr::filter(chr == 5 & position > 74632154-500000 & position < 74657929+500000)
res$target <- "HMGCR"
loci <- rbind(loci, res)

res <- d %>% dplyr::filter(chr == 1 & position > 55505221-500000 & position < 55530525+500000)
res$target <- "PCSK9"
loci <- rbind(loci, res)

res <- d %>% dplyr::filter(chr == 7 & position > 44552134-500000 & position < 44580914+500000)
res$target <- "NPC1L1"
loci <- rbind(loci, res)

res <- d %>% dplyr::filter(chr == 16 & position > 56995762-500000 & position < 57017757+500000)
res$target <- "CETP"
loci <- rbind(loci, res)
d <- loci

# load dosage
snps <- d %>% dplyr::select(chr, position, nea, ea) %>% unique
for (i in 1:nrow(snps)){
    dosage <- tryCatch(
        expr = {
            extract_variant_from_bgen(as.character(snps$chr[i]), as.double(snps$position[i]), snps$nea[i], snps$ea[i])
        },
        error = function(e){ 
            NULL
        }
    )

    if (is.null(dosage)){
        warning(paste0("skipping variant: ", snps$chr[i], ":", snps$position[i]))
        next
    }
    dat <- merge(dat, dosage, "appieu")
}

# drop multiallelics
d <- d %>% filter(key %in% names(dat))
d$trait <- gsub("HDL cholesterol", "hdl_cholesterol.30760.0.0", d$trait)
d$trait <- gsub("LDL cholesterol", "ldl_direct.30780.0.0", d$trait)
d$trait <- gsub("Triglycerides", "triglycerides.30870.0.0", d$trait)

results <- data.frame()
for (i in 1:nrow(d)){
    res <- bp(dat, d$key[i], d$trait[i])
    dat[[paste0(d$trait[i], "_log")]] <- log(dat[[d$trait[i]]])
    res_log <- bp(dat, d$key[i], paste0(d$trait[i], "_log"))
    names(res_log) <- paste0(names(res_log), ".log")
    results <- rbind(results, cbind(res, res_log))
}