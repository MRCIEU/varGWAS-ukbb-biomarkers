library("ieugwasr")
library("dplyr")
library("ggplot2")
library("broom")
source("funs.R")
options(ieugwasr_api="http://64.227.44.193:8006/")
set.seed(123)

f <- "/tmp/tmp.nQ9KZWt0VB/data.33352.csv"
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
#plink \
#--bim /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/exome/released/2020-10-12/data/raw_downloaded/ukb_snp_chr6.bim \
#--fam /projects/MRC-IEU/research/data/ukbiobank/phenotypic/applications/15825/released/2019-05-02/data/data.exome.15825.fam \
#--bed /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/exome/released/2020-10-12/data/raw_downloaded/ukb_cal_chr6_b0_v1.bed \
#--extract glp1r.txt \
#--recode A \
#--out glp1r_exome.txt
dosage <- fread("glp1r_exome.txt.raw")
dosage <- merge(linker, dosage, by.x="app15825", by.y="FID")
dosage <- dosage %>% 
    dplyr::select("appieu", "6:39079018:G:A_G") %>% 
    dplyr::rename(chr6_39046794_G_A="6:39079018:G:A_G") %>%
    dplyr::mutate(chr6_39046794_G_A=dplyr::recode(chr6_39046794_G_A, `2`=0, `0`=2, `1`=1, .default = NA_real_))
dat <- merge(dat, dosage, "appieu", all.x=T)

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

res <- lm(body_mass_index.21001.0.0 ~ chr6_39046794_G_A + age_at_recruitment.21022.0.0 + sex.31.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=dat) %>% tidy %>% dplyr::filter(term == "chr6_39046794_G_A")
res$trait <- "body_mass_index.21001.0.0"
res_mean <- rbind(res_mean, res)
res <- lm(birth_weight_of_first_child.2744.0.0 ~ chr6_39046794_G_A + age_at_recruitment.21022.0.0 + sex.31.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=dat) %>% tidy %>% dplyr::filter(term == "chr6_39046794_G_A")
res$trait <- "birth_weight_of_first_child.2744.0.0"
res_mean <- rbind(res_mean, res)
res <- lm(comparative_body_size_at_age_10.1687.0.0 ~ chr6_39046794_G_A + age_at_recruitment.21022.0.0 + sex.31.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=dat) %>% tidy %>% dplyr::filter(term == "chr6_39046794_G_A")
res$trait <- "comparative_body_size_at_age_10.1687.0.0"
res_mean <- rbind(res_mean, res)

# test for variance effect
res_var <- data.frame()

covar <- c("age_at_recruitment.21022.0.0","sex.31.0.0","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
res <- varGWASR::model(dat %>% dplyr::select(all_of(c("chr6_39016636_C_T", "body_mass_index.21001.0.0", covar))) %>% na.omit, "chr6_39016636_C_T", "body_mass_index.21001.0.0", covar1 = covar, covar2 = covar)
res$trait <- "body_mass_index.21001.0.0"
res$term <- "chr6_39016636_C_T"
res_var <- rbind(res_var, res)

covar <- c("age_at_recruitment.21022.0.0","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
res <- varGWASR::model(dat %>% dplyr::select(all_of(c("chr6_39016636_C_T", "birth_weight_of_first_child.2744.0.0", covar))) %>% na.omit, "chr6_39016636_C_T", "birth_weight_of_first_child.2744.0.0", covar1 = covar, covar2 = covar)
res$trait <- "birth_weight_of_first_child.2744.0.0"
res$term <- "chr6_39016636_C_T"
res_var <- rbind(res_var, res)

covar <- c("age_at_recruitment.21022.0.0","sex.31.0.0","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
res <- varGWASR::model(dat %>% dplyr::select(all_of(c("chr6_39016636_C_T", "comparative_body_size_at_age_10.1687.0.0", covar))) %>% na.omit, "chr6_39016636_C_T", "comparative_body_size_at_age_10.1687.0.0", covar1 = covar, covar2 = covar)
res$trait <- "comparative_body_size_at_age_10.1687.0.0"
res$term <- "chr6_39016636_C_T"
res_var <- rbind(res_var, res)

covar <- c("age_at_recruitment.21022.0.0","sex.31.0.0","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
res <- varGWASR::model(dat %>% dplyr::select(all_of(c("chr6_39046794_G_A", "body_mass_index.21001.0.0", covar))) %>% na.omit, "chr6_39046794_G_A", "body_mass_index.21001.0.0", covar1 = covar, covar2 = covar)
res$trait <- "body_mass_index.21001.0.0"
res$term <- "chr6_39046794_G_A"
res_var <- rbind(res_var, res)

covar <- c("age_at_recruitment.21022.0.0","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
res <- varGWASR::model(dat %>% dplyr::select(all_of(c("chr6_39046794_G_A", "birth_weight_of_first_child.2744.0.0", covar))) %>% na.omit, "chr6_39046794_G_A", "birth_weight_of_first_child.2744.0.0", covar1 = covar, covar2 = covar)
res$trait <- "birth_weight_of_first_child.2744.0.0"
res$term <- "chr6_39046794_G_A"
res_var <- rbind(res_var, res)

covar <- c("age_at_recruitment.21022.0.0","sex.31.0.0","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
res <- varGWASR::model(dat %>% dplyr::select(all_of(c("chr6_39046794_G_A", "comparative_body_size_at_age_10.1687.0.0", covar))) %>% na.omit, "chr6_39046794_G_A", "comparative_body_size_at_age_10.1687.0.0", covar1 = covar, covar2 = covar)
res$trait <- "comparative_body_size_at_age_10.1687.0.0"
res$term <- "chr6_39046794_G_A"
res_var <- rbind(res_var, res)

res_var1 <- res_var %>% dplyr::select(phi_x1, se_x1, trait, term) %>% dplyr::rename(phi=phi_x1, se=se_x1) %>% dplyr::mutate(gt=1)
res_var2 <- res_var %>% dplyr::select(phi_x2, se_x2, trait, term) %>% dplyr::rename(phi=phi_x2, se=se_x2) %>% dplyr::mutate(gt=2)
res_var <- rbind(res_var1, res_var2)
res_var$lci <- res_var$phi - (res_var$se * 1.96)
res_var$uci <- res_var$phi + (res_var$se * 1.96)

# plot
res_var$gt <- as.factor(res_var$gt)
ggplot(res_var, aes(x=gt, y=phi, ymin=lci, ymax=uci)) +
    geom_point() +
    geom_errorbar(width=.05) +
    facet_wrap(trait~term, scales="free_y", ncol = 2, nrow=3) +
    geom_hline(yintercept = c(0), linetype = "dashed", color = "grey") + 
    theme_classic() +
    ylab("Per allele effect on trait (SD) variance (95% CI)") +
    xlab("Copies of GLP1R diabetes decreasing allele")



covar <- c("age_at_recruitment.21022.0.0","sex.31.0.0","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
res <- varGWASR::model(dat %>% dplyr::select(all_of(c("chr6_39046794_G_A", "body_mass_index.21001.0.0_log", covar))) %>% na.omit, "chr6_39046794_G_A", "body_mass_index.21001.0.0_log", covar1 = covar, covar2 = covar)
res$trait <- "body_mass_index.21001.0.0_log"
res$term <- "chr6_39046794_G_A"

covar <- c("age_at_recruitment.21022.0.0","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
res <- varGWASR::model(dat %>% dplyr::select(all_of(c("chr6_39046794_G_A", "birth_weight_of_first_child.2744.0.0_log", covar))) %>% na.omit, "chr6_39046794_G_A", "birth_weight_of_first_child.2744.0.0_log", covar1 = covar, covar2 = covar)
res$trait <- "birth_weight_of_first_child.2744.0.0_log"
res$term <- "chr6_39046794_G_A"

covar <- c("age_at_recruitment.21022.0.0","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
res <- varGWASR::model(dat %>% dplyr::select(all_of(c("chr6_39016636_C_T", "birth_weight_of_first_child.2744.0.0_log", covar))) %>% na.omit, "chr6_39016636_C_T", "birth_weight_of_first_child.2744.0.0_log", covar1 = covar, covar2 = covar)
res$trait <- "birth_weight_of_first_child.2744.0.0_log"
res$term <- "chr6_39016636_C_T"



covar <- c("age_at_recruitment.21022.0.0","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
res <- varGWASR::model(dat %>% dplyr::select(all_of(c("chr6_39046794_G_A", "birth_weight_of_first_child.2744.0.0_log", covar))) %>% na.omit, "chr6_39046794_G_A", "birth_weight_of_first_child.2744.0.0_log", covar1 = covar, covar2 = covar)
res$trait <- "birth_weight_of_first_child.2744.0.0_log"
res$term <- "chr6_39046794_G_A"

covar <- c("age_at_recruitment.21022.0.0","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
res <- varGWASR::model(dat %>% dplyr::select(all_of(c("chr6_39016636_C_T", "birth_weight_of_first_child.2744.0.0_log", covar))) %>% na.omit, "chr6_39016636_C_T", "birth_weight_of_first_child.2744.0.0_log", covar1 = covar, covar2 = covar)
res$trait <- "birth_weight_of_first_child.2744.0.0_log"
res$term <- "chr6_39016636_C_T"