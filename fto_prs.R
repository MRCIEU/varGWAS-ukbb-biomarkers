library("ieugwasr")
library("dplyr")
library("ggplot2")
library("broom")
source("funs.R")
options(ieugwasr_api="http://104.248.164.99/")
set.seed(123)

get_prs <- function(dat, iv){
    # extract dosages and harmonise
    variants <- apply(iv, 1, function(row){
        tryCatch({
            extract_variant_from_bgen(row[['chr']], as.double(row[['position']]), row[['nea']], row[['ea']])
        }, error = function(error_condition) {
            return(NA)
        })
    })
    # drop betas if the variant could not be recovere
    iv <- iv[!is.na(variants),]

    variants <- variants[!is.na(variants)] %>%
        reduce(left_join, by = "appieu")
    row.names(variants) <- variants$appieu
    variants$appieu <- NULL

    # weight each SNP by beta
    w <- as.matrix(variants) %*% diag(iv$beta)

    # sum the weighted alleles into score
    prs <- as.data.frame(rowSums(w))
    names(prs) <- c("PRS")
    prs$appieu <- row.names(prs)

    # merge single snps with pheno
    variants$appieu <- row.names(variants)
    dat <- merge(dat, variants, "appieu")

    # merge prs with pheno
    dat <- merge(dat, prs,"appieu")

    return(dat)
}

f <- "/tmp/tmp.mVvLa5Exif/data.33352.csv"
pheno <- fread(f, select=c(
        "eid",
        "31-0.0",
        "21022-0.0",
        "21001-0.0",
        "2744-0.0",
        "1687-0.0"
    ),
    col.names=c(
        "eid", 
        "sex.31.0.0",
        "age_at_recruitment.21022.0.0",
        "body_mass_index.21001.0.0",
        "birth_weight_of_first_child.2744.0.0",
        "comparative_body_size_at_age_10.1687.0.0")
)

pheno <- pheno %>% dplyr::mutate_at(c('birth_weight_of_first_child.2744.0.0'), na_if, -1)
pheno <- pheno %>% dplyr::mutate_at(c('birth_weight_of_first_child.2744.0.0'), na_if, -2)
pheno <- pheno %>% dplyr::mutate_at(c('birth_weight_of_first_child.2744.0.0'), na_if, -3)
pheno <- pheno %>% dplyr::mutate_at(c('comparative_body_size_at_age_10.1687.0.0'), na_if, -1)
pheno <- pheno %>% dplyr::mutate_at(c('comparative_body_size_at_age_10.1687.0.0'), na_if, -3)
pheno <- pheno %>% dplyr::mutate(comparative_body_size_at_age_10.1687.0.0=dplyr::recode(comparative_body_size_at_age_10.1687.0.0, `3`=0, `1`=-1, `2`=1, .default = NA_real_))
pheno$birth_weight_of_first_child.2744.0.0 <- as.double(pheno$birth_weight_of_first_child.2744.0.0)

# load linker
linker <- get_filtered_linker(drop_standard_excl=TRUE, drop_non_white_british=TRUE, drop_related=TRUE, application="15825")

# load covariates
pc <- get_genetic_principal_components()

# merge data
dat <- merge(linker, pc, "appieu")
dat <- merge(dat, pheno, by.x="app15825", by.y="eid")

# SD scale outcomes
dat$body_mass_index.21001.0.0 <- dat$body_mass_index.21001.0.0 / sd(dat$body_mass_index.21001.0.0, na.rm=T)
dat$birth_weight_of_first_child.2744.0.0 <- dat$birth_weight_of_first_child.2744.0.0 / sd(dat$birth_weight_of_first_child.2744.0.0, na.rm=T)
dat$comparative_body_size_at_age_10.1687.0.0 <- dat$comparative_body_size_at_age_10.1687.0.0 / sd(dat$comparative_body_size_at_age_10.1687.0.0, na.rm=T)

# extract array SNP
dosage <- extract_variant_from_bgen("16", 53803574, "T", "A")
dat <- merge(dat, dosage, "appieu")

# prepare BMI-PRS lacking GLP1R locus
bmi_iv <- tophits("ieu-b-40")
bmi_dat <- get_prs(dat, bmi_iv %>% dplyr::filter(chr != "16"))

# test for association of FTO with BMI along quntiles of the BMI-PRS
