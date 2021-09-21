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

option_list = list(
  make_option(c("-t", "--trait"), type="character", default=NULL, help="Name of trait", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

message(paste0("trait ", opt$trait))

# read in extracted phenotypes
pheno <- fread(paste0("data/", opt$trait, ".txt"))

# read top gxg effects on additive scale and filter out effects without support on multiplicative scale
d <- fread("data/gxg.txt")
d <- d %>% filter(trait != "body_mass_index.21001.0.0")
d <- d %>% filter(p.value < 5e-8)
d <- cbind(d, as.data.frame(str_split(d$term, ":", simplify=T), stringsAsFactors=F))
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

# test for effect of SNP stratified by modifier
results <- data.frame()
for (i in 1:nrow(d)){
  if (d$trait[i] != opt$trait){
    next
  }
  # test SNP effect by stratified modifier
  f <- as.formula(paste0(opt$trait, " ~ ", d$V1[i], " + age_at_recruitment.21022.0.0 + sex.31.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))
  k <- d$V2[i]
  pheno$mod_gt <- round(pheno[[k]])
  pheno0 <- pheno %>% dplyr::filter(mod_gt == 0)
  pheno1 <- pheno %>% dplyr::filter(mod_gt == 1)
  pheno2 <- pheno %>% dplyr::filter(mod_gt == 2)

  mod <- lm(f, data=pheno0)
  t <- coeftest(mod, vcov = vcovHC(mod, type = "HC0")) %>% tidy
  t <- t %>% dplyr::filter(d$V1[i] == t$term)
  t$mod <- 0
  t$mod_snp <- k
  t$trait <- opt$trait
  results <- rbind(results, t)

  mod <- lm(f, data=pheno1)
  t <- coeftest(mod, vcov = vcovHC(mod, type = "HC0")) %>% tidy
  t <- t %>% dplyr::filter(d$V1[i] == t$term)
  t$mod <- 1
  t$mod_snp <- k
  t$trait <- opt$trait
  results <- rbind(results, t)

  mod <- lm(f, data=pheno2)
  t <- coeftest(mod, vcov = vcovHC(mod, type = "HC0")) %>% tidy
  t <- t %>% dplyr::filter(d$V1[i] == t$term)
  t$mod <- 2
  t$mod_snp <- k
  t$trait <- opt$trait
  results <- rbind(results, t)

}

# save
write.table(results, sep="\t", quote=F, row.names=F, file=paste0("data/", opt$trait, ".gxg-qual.txt"))