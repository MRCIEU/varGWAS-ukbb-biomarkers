load("data/pheno.RData")
library('optparse')
library('data.table')
library('dplyr')
library('broom')
library('ieugwasr')
library("plm")
library('sandwich')
library('lmtest')
source("funs.R")
set.seed(1234)

# Adapted from: https://github.com/nmdavies/within_family_mr/blob/master/HUNT/reg_2_option_c_individual_snps-v5.R
related_plm <- function(df, out, snp1, snp2){   
    #Run plm model with family fixed effects + robust standard errors
    f <- as.formula(paste(out, "~ ",snp1,"*", snp2," + age_at_recruitment.21022.0.0 + sex.31.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))
    estimates<-plm(f, data = df, index = "famid", model = "within", inst.method = "bvk")
    estimates_robust<-coeftest(estimates,vcov=vcovHC(estimates,type="HC0",cluster="group"))
    
    mask <- grepl(":" ,row.names(estimates_robust))
    beta<-estimates_robust[mask,1]
    se<-estimates_robust[mask,2]
    pvalue<-estimates_robust[mask,4]
    sample_size<-nrow(model.frame(estimates))

    return(data.frame(beta, se, pvalue, sample_size, term=row.names(estimates_robust)[mask]))
}

get_siblings <- function(){
    dat <- fread("/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/derived/related/siblings/Siblings.01.fam")
    names(dat)[1]<-'famid'
    names(dat)[2]<-'appieu'
    return(dat[,1:2])
}

option_list = list(
  make_option(c("-t", "--trait"), type="character", default=NULL, help="Name of trait", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

message(paste0("trait ", opt$trait))

# load linker
linker <- get_filtered_linker(drop_standard_excl=TRUE, drop_non_white_british=FALSE, drop_related=FALSE, application="15825")

# load covariates
covariates <- get_covariates()
covariates$chip <- as.numeric(as.factor(covariates$chip)) - 1
pc <- get_genetic_principal_components()

# merge data
sibs <- get_siblings()
linker <- merge(linker, sibs, "appieu")
dat <- merge(linker, pheno, by.x="app15825", by.y="eid")
dat <- merge(dat, covariates, "appieu")
dat <- merge(dat, pc, "appieu")

# select fields for GWAS
dat <- dat[,c("appieu", "sex.31.0.0", "age_at_recruitment.21022.0.0", "famid", "chip", opt$trait, paste0("PC", seq(1, 20))), with=F]

# drop missing values
dat <- dat[complete.cases(dat), ]

# SD scale
dat[[opt$trait]] <- dat[[opt$trait]] / sd(dat[[opt$trait]])
dat$age_at_recruitment.21022.0.0 <- dat$age_at_recruitment.21022.0.0 / sd(dat$age_at_recruitment.21022.0.0)
pheno <- dat
pheno$famid <- as.factor(pheno$famid)

# read in vGWAS
snps <- get_variants(opt$trait)

# filter on P value
vqtls <- snps %>% filter(phi_p < 5e-5) %>% select("rsid", "phi_p") %>% rename(pval = phi_p)

# clump records
vqtls <- ld_clump(vqtls)
snps <- snps[snps$rsid %in% vqtls$rsid]

# add key
snps$key <- paste0("chr", snps$chr, "_", snps$pos, "_", snps$oa, "_", snps$ea)

# load dosages
for (i in 1:nrow(snps)){
    dosage <- extract_variant_from_bgen(as.character(snps$chr[i]), as.double(snps$pos[i]), snps$oa[i], snps$ea[i])
    pheno <- merge(pheno, dosage, "appieu")
}

# select vQTLs
vqtls <- grep("^chr", names(pheno), value=T)

# drop missing values
pheno <- pheno[complete.cases(pheno), ]

# drop singletons
pheno <- pheno %>%
  group_by(famid) %>%
  filter(n() > 1)

# test for interaction between each snp
results <- data.frame()
for (i in 1:length(vqtls)){
  for (j in 1:length(vqtls)){

    # skip GxG on same chromosome within 10Mb
    i_chr <- snps %>% filter(key == vqtls[i]) %>% pull("chr")
    i_pos <- snps %>% filter(key == vqtls[i]) %>% pull("pos")
    j_chr <- snps %>% filter(key == vqtls[j]) %>% pull("chr")
    j_pos <- snps %>% filter(key == vqtls[j]) %>% pull("pos")

    if (i_chr == j_chr){
      if (abs(i_pos - j_pos) < 10000000){
        message("Skipping test (<10Mb) for: ", vqtls[i], " ", vqtls[j])
        next
      }
    }

    if (paste0(vqtls[j], ":" ,vqtls[i]) %in% results$term){
        message("Skipping test (already done) for: ", vqtls[i], " ", vqtls[j])
        next
    }

    # test GxG
    message("Testing GxG for: ", vqtls[i], " ", vqtls[j])
    fit <- related_plm(pheno, opt$trait, vqtls[i], vqtls[j])

    # store results
    results <- rbind(results, fit)
  }
}

# save
write.table(results, sep="\t", quote=F, row.names=F, file=paste0("data/", opt$trait, ".gxg-wf.txt"))