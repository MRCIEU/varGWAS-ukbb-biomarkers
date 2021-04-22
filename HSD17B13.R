library('data.table')
library('dplyr')
library('ieugwasr')
library('broom')
library('TwoSampleMR')
source("funs.R")
set.seed(124)

### MR effect of ALT on HSD17B13 expression ###
# Get instruments
exposure_dat <- extract_instruments("ukb-d-30620_irnt")

# Get effects of instruments on outcome
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="eqtl-a-ENSG00000170509")

# Harmonise the exposure and outcome data
dat <- harmonise_data(exposure_dat, outcome_dat)

# Perform MR
mr_alt_hsd <- mr(dat)

### MR effect of HSD17B13 expression on ALT ###
# Get instruments
exposure_dat <- extract_instruments("eqtl-a-ENSG00000170509")

# Get effects of instruments on outcome
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ukb-d-30620_irnt")

# Harmonise the exposure and outcome data
dat <- harmonise_data(exposure_dat, outcome_dat)

# Perform MR
mr_hsd_alt <- mr(dat)

### LoF effect of HSD17B13 expression on ALT ###

# read in extracted phenotypes
pheno <- fread("data/alanine_aminotransferase.30620.txt")

# read in LoF calls
calls <- fread("data/ukbb.HSD17B13.exome.txt.raw")

# link to IEU
linker <- get_filtered_linker(drop_standard_excl=FALSE, drop_non_white_british=FALSE, drop_related=FALSE, app="15825")
calls <- merge(calls, linker, by.x="FID", by.y="app15825")

# load variant annotations
anno <- fread("data/ukbb.HSD17B13.exome.vep.txt", skip=42, check=T)
anno$stop <- grepl("frameshift_variant|stop_gained", anno$Consequence)

# select variants
high_impact <- anno %>%
    filter(Feature == "ENST00000328546") %>%
    filter(IMPACT == "HIGH") %>%
    select("X.Uploaded_variation") %>%
    distinct() %>%
    pull()

stop <- anno %>%
    filter(Feature == "ENST00000328546") %>%
    filter(stop) %>%
    select("X.Uploaded_variation") %>%
    distinct() %>%
    pull()

high_impact_missense <- anno %>%
    filter(Feature == "ENST00000328546") %>%
    filter(IMPACT == "HIGH" | (grepl("missense", Consequence) & SIFT == "deleterious" & PolyPhen == "probably_damaging")) %>%
    select("X.Uploaded_variation") %>%
    distinct() %>%
    pull()

# select high impact calls from total callset
high_impact_calls <- calls[,grepl(paste(high_impact, collapse="|"), names(calls)), with=F]
high_impact_calls <- data.frame(appieu=calls$appieu, n_high_impact=rowSums(high_impact_calls, na.rm=T))
high_impact_calls$has_high_impact <- high_impact_calls$n_high_impact != 0

stop_calls <- calls[,grepl(paste(stop, collapse="|"), names(calls)), with=F]
stop_calls <- data.frame(appieu=calls$appieu, n_stop=rowSums(stop_calls, na.rm=T))
stop_calls$has_stop <- stop_calls$n_stop != 0

high_impact_missense_calls <- calls[,grepl(paste(high_impact_missense, collapse="|"), names(calls)), with=F]
high_impact_missense_calls <- data.frame(appieu=calls$appieu, n_high_impact_missense=rowSums(high_impact_missense_calls, na.rm=T))
high_impact_missense_calls$has_high_impact_missense <- high_impact_missense_calls$n_high_impact_missense != 0

# merge
pheno <- merge(pheno, high_impact_calls, "appieu")
pheno <- merge(pheno, high_impact_missense_calls, "appieu")
pheno <- merge(pheno, stop_calls, "appieu")

# LoF effect on ALT mean
fits <- lm(alanine_aminotransferase.30620 ~ has_stop + age_at_recruitment.21022.0.0 + sex.31.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=pheno)
fitm <- lm(alanine_aminotransferase.30620 ~ has_high_impact_missense + age_at_recruitment.21022.0.0 + sex.31.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=pheno)
fith <- lm(alanine_aminotransferase.30620 ~ has_high_impact + age_at_recruitment.21022.0.0 + sex.31.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=pheno)

# missense effect on ALT variance
fit1 <- lm(alanine_aminotransferase.30620 ~ has_high_impact_missense + age_at_recruitment.21022.0.0 + sex.31.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=pheno)
pheno$d <- resid(fit1)^2
fit2 <- lm(d ~ has_high_impact_missense, data=pheno)