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

# ALT
snps <- ieugwasr::tophits("ukb-d-30620_irnt")
for (i in 1:nrow(snps)){
    dosage <- tryCatch(
        expr = {
            extract_variant_from_bgen(as.character(snps$chr[i]), as.double(snps$position[i]), snps$nea[i], snps$ea[i])
        },
        error = function(e){ 
            message("skipping:", paste0(as.character(snps$chr[i]), as.double(snps$position[i]), snps$nea[i], snps$ea[i]))
            NULL
        }
    )
    if (!is.null(dosage)){
        dat <- merge(dat, dosage, "appieu")
    }
}
# test for subgroup effect
f <- paste0("alanine_aminotransferase.30620.0.0 ~ chr9_132566666_G_A + age_at_recruitment.21022.0.0 + sex.31.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + ", paste0(grep("chr", names(dat), value=T), collapse="+"))
as.formula(f)
fit <- lm(f, data=dat) %>% tidy
fit1 <- lm(f, data=dat %>% dplyr::filter(body_mass_index.21001.0.0_b)) %>% tidy
fit0 <- lm(f, data=dat %>% dplyr::filter(!body_mass_index.21001.0.0_b)) %>% tidy
fit %>% dplyr::filter(term == "chr9_132566666_G_A")
fit0 %>% dplyr::filter(term == "chr9_132566666_G_A")
fit1 %>% dplyr::filter(term == "chr9_132566666_G_A")

# TG
snps <- ieugwasr::tophits("ukb-d-30870_irnt")
for (i in 1:nrow(snps)){
    dosage <- tryCatch(
        expr = {
            extract_variant_from_bgen(as.character(snps$chr[i]), as.double(snps$position[i]), snps$nea[i], snps$ea[i])
        },
        error = function(e){ 
            message("skipping:", paste0(as.character(snps$chr[i]), as.double(snps$position[i]), snps$nea[i], snps$ea[i]))
            NULL
        }
    )
    if (!is.null(dosage)){
        dat <- merge(dat, dosage, "appieu")
    }
}
dosage <- extract_variant_from_bgen("22", 44324727, "C", "G")
dat <- merge(dat, dosage, "appieu")
# test for subgroup effect
f <- paste0("triglycerides.30870.0.0 ~ chr22_44324727_C_G + age_at_recruitment.21022.0.0 + sex.31.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + ", paste0(grep("chr", names(dat), value=T), collapse="+"))
as.formula(f)
fit <- lm(f, data=dat) %>% tidy
fit1 <- lm(f, data=dat %>% dplyr::filter(body_mass_index.21001.0.0_b)) %>% tidy
fit0 <- lm(f, data=dat %>% dplyr::filter(!body_mass_index.21001.0.0_b)) %>% tidy
fit %>% dplyr::filter(term == "chr22_44324727_C_G")
fit0 %>% dplyr::filter(term == "chr22_44324727_C_G")
fit1 %>% dplyr::filter(term == "chr22_44324727_C_G")

# Urate
snps <- ieugwasr::tophits("ukb-d-30880_irnt")
for (i in 1:nrow(snps)){
    dosage <- tryCatch(
        expr = {
            extract_variant_from_bgen(as.character(snps$chr[i]), as.double(snps$position[i]), snps$nea[i], snps$ea[i])
        },
        error = function(e){ 
            message("skipping:", paste0(as.character(snps$chr[i]), as.double(snps$position[i]), snps$nea[i], snps$ea[i]))
            NULL
        }
    )
    if (!is.null(dosage)){
        dat <- merge(dat, dosage, "appieu")
    }
}
dosage <- extract_variant_from_bgen("4", 10402838, "T", "C")
dat <- merge(dat, dosage, "appieu")
# test for subgroup effect
f <- paste0("urate.30880.0.0 ~ chr4_10402838_T_C + age_at_recruitment.21022.0.0 +  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + ", paste0(grep("chr", names(dat), value=T), collapse="+"))
as.formula(f)
fit <- lm(f, data=dat) %>% tidy
fit1 <- lm(f, data=dat %>% dplyr::filter(sex.31.0.0_b)) %>% tidy
fit0 <- lm(f, data=dat %>% dplyr::filter(!sex.31.0.0_b)) %>% tidy
fit %>% dplyr::filter(term == "chr4_10402838_T_C")
fit0 %>% dplyr::filter(term == "chr4_10402838_T_C")
fit1 %>% dplyr::filter(term == "chr4_10402838_T_C")