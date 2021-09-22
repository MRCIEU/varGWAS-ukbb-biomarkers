library('optparse')
library('data.table')
library('dplyr')
library('broom')
library('ieugwasr')
library("lmtest")
library("sandwich")
source("funs.R")
set.seed(1234)
options(ieugwasr_api="http://64.227.44.193:8006/")

# get liver disease outcomes
load("data/pheno.RData")
linker <- get_filtered_linker(drop_standard_excl=TRUE, drop_non_white_british=TRUE, drop_related=TRUE, application="15825")
dat <- merge(linker, pheno, by.x="app15825", by.y="eid")
dat <- dat %>% select(appieu, liver_disease, alcoholic_liver_disease, fibrosis_liver_disease, fatty_liver_disease)

# read in extracted phenotypes & merge with liver disease
pheno <- fread(paste0("data/alanine_aminotransferase.30620.0.0.txt"))
pheno <- merge(pheno, dat, "appieu")

# read in ALT effects PNPLA3 x HSD17B13
d <- fread("data/gxg.txt")
d <- d %>% filter(trait != "body_mass_index.21001.0.0")
d <- d %>% filter(p.value < 5e-8)
d <- cbind(d, as.data.frame(str_split(d$term, ":", simplify=T), stringsAsFactors=F))
d$V3 <- d$V1
d$V1 <- d$V2
d$V2 <- d$V3
d$V3 <- NULL
v1 <- as.data.frame(str_split(d$V1, "_", simplify=T), stringsAsFactors=F)
names(v1) <- c("chr.1", "pos.1", "oa.1", "ea.1")
v2 <- as.data.frame(str_split(d$V2, "_", simplify=T), stringsAsFactors=F)
names(v2) <- c("chr.2", "pos.2", "oa.2", "ea.2")
d <- cbind(d, v1, v2)
d$chr.1 <- gsub("chr", "", d$chr.1)
d$chr.2 <- gsub("chr", "", d$chr.2)

# load dosages
snps <- d %>% select(chr.1, pos.1, oa.1, ea.1) %>% rename(chr="chr.1", pos="pos.1", oa="oa.1", ea="ea.1")
snps <- rbind(snps, d %>% select(chr.2, pos.2, oa.2, ea.2) %>% rename(chr="chr.2", pos="pos.2", oa="oa.2", ea="ea.2"))
snps <- unique(snps)
snps$key <- paste0("chr", snps$chr, "_", snps$pos, "_", snps$oa, "_", snps$ea)

for (i in 1:nrow(snps)){
    dosage <- extract_variant_from_bgen(as.character(snps$chr[i]), as.double(snps$pos[i]), snps$oa[i], snps$ea[i])
    pheno <- merge(pheno, dosage, "appieu")
}

stopifnot(nrow(d) == 1)

#stratify
k <- "chr4_88212722_A_G" # HSDB1317
pheno$mod_gt <- round(pheno[[k]])
pheno0 <- pheno %>% dplyr::filter(mod_gt == 0)
pheno1 <- pheno %>% dplyr::filter(mod_gt == 1)
pheno2 <- pheno %>% dplyr::filter(mod_gt == 2)

# test for effect of SNP stratified by modifier
results <- data.frame()

# test SNP effect on ALT by stratified modifier
f <- as.formula(paste0("alanine_aminotransferase.30620.0.0", " ~ chr22_44324730_T_C + age_at_recruitment.21022.0.0 + sex.31.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))
mod <- lm(f, data=pheno0)
t <- coeftest(mod, vcov = vcovHC(mod, type = "HC0")) %>% tidy
t <- t[2,]
t$mod <- 0
t$mod_snp <- k
t$trait <- "alanine_aminotransferase.30620.0.0"
results <- rbind(results, t)

mod <- lm(f, data=pheno1)
t <- coeftest(mod, vcov = vcovHC(mod, type = "HC0")) %>% tidy
t <- t[2,]
t$mod <- 1
t$mod_snp <- k
t$trait <- "alanine_aminotransferase.30620.0.0"
results <- rbind(results, t)

mod <- lm(f, data=pheno2)
t <- coeftest(mod, vcov = vcovHC(mod, type = "HC0")) %>% tidy
t <- t[2,]
t$mod <- 2
t$mod_snp <- k
t$trait <- "alanine_aminotransferase.30620.0.0"
results <- rbind(results, t)

# test SNP effect on liver disease by stratified modifier
for (outcome in c("alcoholic_liver_disease", "fibrosis_liver_disease", "fatty_liver_disease")){
  f <- as.formula(paste0(outcome, " ~ chr22_44324730_T_C + age_at_recruitment.21022.0.0 + sex.31.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))
  mod <- glm(f, data=pheno0, family="binomial")
  t <- tidy(mod)[2,]
  t$mod <- 0
  t$mod_snp <- k
  t$trait <- outcome
  results <- rbind(results, t)

  mod <- glm(f, data=pheno1, family="binomial")
  t <- tidy(mod)[2,]
  t$mod <- 1
  t$mod_snp <- k
  t$trait <- outcome
  results <- rbind(results, t)

  mod <- glm(f, data=pheno2, family="binomial")
  t <- tidy(mod)[2,]
  t$mod <- 2
  t$mod_snp <- k
  t$trait <- outcome
  results <- rbind(results, t)
}

# save
write.table(results, sep="\t", quote=F, row.names=F, file=paste0("data/alanine_aminotransferase.30620.0.0.gxg-qual.txt"))