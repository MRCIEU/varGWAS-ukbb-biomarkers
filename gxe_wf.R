load("data/pheno.RData")
library('optparse')
library('data.table')
library('dplyr')
library('broom')
library('ieugwasr')
library("plm")
library('sandwich')
library('lmtest')
library("stringr")
library('forestplot')
source("funs.R")
set.seed(1234)
options(ieugwasr_api="http://64.227.44.193:8006/")

option_list = list(
  make_option(c("-t", "--trait"), type="character", default=NULL, help="Name of trait", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

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
dat <- dat[,c("appieu", "famid", "chip", opt$trait, env_exp, paste0("PC", seq(1, 20))), with=F]

# drop missing values
dat <- dat[complete.cases(dat), ]

# SD scale
dat[[opt$trait]] <- dat[[opt$trait]] / sd(dat[[opt$trait]])
# SD scale env
for (e in env_exp){
  dat[[e]] <- dat[[e]] / sd(dat[[e]], na.rm=T)
}
pheno <- dat
pheno$famid <- as.factor(pheno$famid)

# read in GxE associations
snps <- fread(paste0("data/", opt$trait, ".gxe.txt"))

# take top GxE hits for replication
snps <- snps %>% filter(p.value < 5e-8)

# split term
snps <- cbind(snps, str_split(snps$term, ":", simplify=T), stringsAsFactors=F)
usnps <- as.data.frame(str_split(snps$V2, "_", simplify=T), stringsAsFactors=F)
usnps$V1 <- gsub("chr", "", usnps$V1)
usnps$V2 <- as.numeric(usnps$V2)
usnps <- unique(usnps)

# load dosages
for (i in 1:nrow(usnps)){
    dosage <- extract_variant_from_bgen(usnps$V1[i], usnps$V2[i], usnps$V3[i], usnps$V4[i])
    pheno <- merge(pheno, dosage, "appieu")
}

# drop missing values
pheno <- pheno[complete.cases(pheno), ]

# drop singletons
pheno <- pheno %>%
  group_by(famid) %>%
  filter(n() > 1)

# test for interaction between each snp
results <- data.frame()
for (i in 1:length(snps$term)){
  # test GxE
  pair <- str_split(snps$term[i], ":", simplify=T)
  message("Testing GxE for: ", pair[1], " ", pair[2])
  fit <- related_plm(pheno, opt$trait, pair[1], pair[2])

  # store results
  results <- rbind(results, fit)
}

# save
write.table(results, sep="\t", quote=F, row.names=F, file=paste0("data/", opt$trait, ".gxe-wf.txt"))