load("data/pheno.RData")
library("data.table")
library("dplyr")
library("stringr")
library("broom")
library('optparse')
source("funs.R")
set.seed(123)

option_list = list(
  make_option(c("-t", "--trait"), type="character", default=NULL, help="Name of trait", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# load linker
linker <- get_filtered_linker(drop_standard_excl=TRUE, drop_non_white_british=TRUE, drop_related=TRUE, application="15825")

# load covariates
covariates <- get_covariates()
covariates$chip <- as.numeric(as.factor(covariates$chip)) - 1
pc <- get_genetic_principal_components()

# merge data
dat <- merge(linker, pheno, by.x="app15825", by.y="eid")
dat <- merge(dat, covariates, "appieu")
dat <- merge(dat, pc, "appieu")

# read in GxG associations
snps <- fread(paste0("data/", opt$trait, ".gxg.txt"))

# filter on P value
snps <- snps %>% filter(p.value < 5e-8)

# split term
snps <- cbind(snps, str_split(snps$term, ":", simplify=T), stringsAsFactors=F)
v1 <- as.data.frame(str_split(snps$V1, "_", simplify=T), stringsAsFactors=F)
names(v1) <- c("chr.1", "pos.1", "oa.1", "ea.1")
v2 <- as.data.frame(str_split(snps$V2, "_", simplify=T), stringsAsFactors=F)
names(v2) <- c("chr.2", "pos.2", "oa.2", "ea.2")
snps <- cbind(snps, v1, v2, stringsAsFactors=F)
snps$chr.1 <- gsub("chr", "", snps$chr.1)
snps$chr.2 <- gsub("chr", "", snps$chr.2)
snps$pos.1 <- as.numeric(snps$pos.1)
snps$pos.2 <- as.numeric(snps$pos.2)

# test for effect on disease outcomes
results <- data.frame()
for (i in 1:nrow(snps)){
  snp1 <- extract_variant_from_bgen(snps$chr.1[i], snps$pos.1[i], snps$oa.1[i], snps$ea.1[i])
  snp2 <- extract_variant_from_bgen(snps$chr.2[i], snps$pos.2[i], snps$oa.2[i], snps$ea.2[i])
  temp <- merge(dat, snp1, "appieu")
  temp <- merge(temp, snp2, "appieu")
  for (out in c("liver_disease", "fatty_liver_disease", "fibrosis_liver_disease", "alcoholic_liver_disease", "CKD", "gout", "T2DM", "heart_attack.6150", "stroke.6150")){
    f <- as.formula(paste0(out, " ~ chr", snps$chr.1[i], "_", snps$pos.1[i], "_",  snps$oa.1[i], "_", snps$ea.1[i], "* chr", snps$chr.2[i], "_", snps$pos.2[i], "_",  snps$oa.2[i], "_", snps$ea.2[i], " + age_at_recruitment.21022.0.0 + sex.31.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))
    fit <- glm(f, data=temp)
    fit <- tidy(fit) %>% filter(grepl(":", term))
    fit$out <- out
    results <- rbind(results, fit)
  }
}

write.table(results, sep="\t", quote=F, row.names=F, file=paste0("data/", opt$trait, ".gxg-disease.txt"))