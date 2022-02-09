library("ieugwasr")
library("dplyr")
library("ggplot2")
library("broom")
source("funs.R")
options(ieugwasr_api="http://64.227.44.193:8006/")
set.seed(123)

f <- "/tmp/tmp.ALUSQjx5OJ/data.33352.csv"
pheno <- fread(f, select=c(
        "eid",
        "31-0.0",
        "21022-0.0",
        "21001-0.0",
        "4080-0.0"
    ),
    col.names=c(
        "eid", 
        "sex.31.0.0",
        "age_at_recruitment.21022.0.0",
        "body_mass_index.21001.0.0",
        "systolic_blood_pressure_automated_reading.4080.0.0")
)

# load linker
linker <- get_filtered_linker(drop_standard_excl=TRUE, drop_non_white_british=TRUE, drop_related=TRUE, application="15825")

# load covariates
pc <- get_genetic_principal_components()

# merge data
dat <- merge(linker, pc, "appieu")
dat <- merge(dat, pheno, by.x="app15825", by.y="eid")

# SD scale outcomes
dat$body_mass_index.21001.0.0 <- dat$body_mass_index.21001.0.0 / sd(dat$body_mass_index.21001.0.0, na.rm=T)
dat$systolic_blood_pressure_automated_reading.4080.0.0 <- dat$systolic_blood_pressure_automated_reading.4080.0.0 / sd(dat$systolic_blood_pressure_automated_reading.4080.0.0, na.rm=T)

# select SNPs and extract from UKBB
bmi_iv <- ieugwasr::tophits("ieu-a-835")
for (i in 1:nrow(bmi_iv)){
    dosage <- extract_variant_from_bgen(bmi_iv$chr[i], as.double(bmi_iv$position[i]), bmi_iv$nea[i], bmi_iv$ea[i])
    dat <- merge(dat, dosage, "appieu")
}

# test for SNP-variance effects
covar <- c("age_at_recruitment.21022.0.0","sex.31.0.0","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
res_var <- data.frame()
for (i in 1:nrow(bmi_iv)){
    snp <- paste0("chr", bmi_iv$chr[i], "_", as.double(bmi_iv$position[i]), "_", bmi_iv$nea[i], "_", bmi_iv$ea[i])
    
    # BMI
    res <- varGWASR::model(dat %>% dplyr::select(all_of(c(snp, "body_mass_index.21001.0.0", covar))) %>% na.omit, snp, "body_mass_index.21001.0.0", covar1 = covar, covar2 = covar)
    res$trait <- "body_mass_index.21001.0.0"
    res$term <- bmi_iv$rsid[i]
    res$ea <- bmi_iv$ea[i]
    res$nea <- bmi_iv$nea[i]
    res_var <- rbind(res_var, res)

    # SBP
    res <- varGWASR::model(dat %>% dplyr::select(all_of(c(snp, "systolic_blood_pressure_automated_reading.4080.0.0", covar))) %>% na.omit, snp, "systolic_blood_pressure_automated_reading.4080.0.0", covar1 = covar, covar2 = covar)
    res$trait <- "systolic_blood_pressure_automated_reading.4080.0.0"
    res$term <- bmi_iv$rsid[i]
    res$ea <- bmi_iv$ea[i]
    res$nea <- bmi_iv$nea[i]
    res_var <- rbind(res_var, res)
}

# save file
write.table(res_var, sep="\t", quote=F, row.names=F, file="bmi_iv_var.txt")