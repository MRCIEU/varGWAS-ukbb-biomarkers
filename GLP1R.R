library("ieugwasr")
library("dplyr")
library("broom")
source("funs.R")
options(ieugwasr_api="http://64.227.44.193:8006/")
set.seed(123)

f <- "/tmp/tmp.GU8Wh7SknM/data.33352.csv"
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
dosage <- extract_variant_from_bgen("6", 39016636, "C", "T")
dat <- merge(dat, dosage, "appieu")

# extract exome SNP
# TODO
#system("
#    plink \
#    --bim /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/exome/released/2020-10-12/data/raw_downloaded/ukb_snp_chr6.bim \
#    --fam /projects/MRC-IEU/research/data/ukbiobank/genetic/variants/exome/raw/exome_download/data.exome.ukbapp.15825.fam \
#    --bed /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/exome/released/2020-10-12/data/raw_downloaded/ukb_cal_chr6_b0_v1.bed \
#    --extract <(grep -v '#' ~/ukbb.lct.exome.vep.txt | cut -s -f1 | sort -u) \
#    --recode A \
#    --out ~/ukbb.lct.exome.txt
#")

# test for mean effect
res_mean <- data.frame()

res <- lm(body_mass_index.21001.0.0 ~ chr6_39016636_C_T + age_at_recruitment.21022.0.0 + sex.31.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=dat) %>% tidy %>% dplyr::filter(term == "chr6_39016636_C_T")
res$trait <- "body_mass_index.21001.0.0"
res_mean <- rbind(res_mean, res)
res <- lm(birth_weight_of_first_child.2744.0.0 ~ chr6_39016636_C_T + age_at_recruitment.21022.0.0 + sex.31.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=dat) %>% tidy %>% dplyr::filter(term == "chr6_39016636_C_T")
res$trait <- "birth_weight_of_first_child.2744.0.0"
res_mean <- rbind(res_mean, res)
res <- lm(comparative_body_size_at_age_10.1687.0.0 ~ chr6_39016636_C_T + age_at_recruitment.21022.0.0 + sex.31.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=dat) %>% tidy %>% dplyr::filter(term == "chr6_39016636_C_T")
res$trait <- "comparative_body_size_at_age_10.1687.0.0"
res_mean <- rbind(res_mean, res)

# test for variance effect
res_var <- data.frame()

covar <- c("age_at_recruitment.21022.0.0","sex.31.0.0","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
res <- varGWASR::model(dat %>% dplyr::select(all_of(c("chr6_39016636_C_T", "body_mass_index.21001.0.0", covar))) %>% na.omit, "chr6_39016636_C_T", "body_mass_index.21001.0.0", covar1 = covar, covar2 = covar)
res$trait <- "body_mass_index.21001.0.0"
res_var <- rbind(res_var, res)

covar <- c("age_at_recruitment.21022.0.0","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
res <- varGWASR::model(dat %>% dplyr::select(all_of(c("chr6_39016636_C_T", "birth_weight_of_first_child.2744.0.0", covar))) %>% na.omit, "chr6_39016636_C_T", "birth_weight_of_first_child.2744.0.0", covar1 = covar, covar2 = covar)
res$trait <- "birth_weight_of_first_child.2744.0.0"
res_var <- rbind(res_var, res)

covar <- c("age_at_recruitment.21022.0.0","sex.31.0.0","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
res <- varGWASR::model(dat %>% dplyr::select(all_of(c("chr6_39016636_C_T", "comparative_body_size_at_age_10.1687.0.0", covar))) %>% na.omit, "chr6_39016636_C_T", "comparative_body_size_at_age_10.1687.0.0", covar1 = covar, covar2 = covar)
res$trait <- "comparative_body_size_at_age_10.1687.0.0"
res_var <- rbind(res_var, res)