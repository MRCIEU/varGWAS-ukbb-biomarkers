library('optparse')
library('data.table')
library('broom')
source("funs.R")
set.seed(1234)

option_list = list(
  make_option(c("-p", "--pheno_file"), type="character", default=NULL, help="UKBB pheno CSV", metavar="character"),
  make_option(c("-t", "--trait"), type="character", default=NULL, help="Name of trait", metavar="character"),
  make_option(c("-s", "--snp_file"), type="character", default=NULL, help="SNP file", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, help="Output file", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

opt <- data.frame(p="data/urate.30880.0.0.txt", t="urate.30880.0.0", s="data/urate.30880.0.0.clump.txt", stringsAsFactors=F)

# read in extracted phenotypes
pheno <- fread(opt$p)

# read in snp list
snps <- fread(opt$s)

# load dosages
for (i in 1:nrow(snps)){
    dosage <- extract_variant_from_bgen(as.character(snps$chr[i]), as.double(snps$pos[i]), snps$oa[i], snps$ea[i])
    pheno <- merge(pheno, dosage, "appieu")
}

# select vQTLs
vqtls <- grep("^chr", names(pheno), value=T)

# test for interaction between each snp
results <- data.frame()
for (i in 1:length(vqtls)){
  f <- as.formula(paste0(opt$t, " ~ age_at_recruitment.21022.0.0 + sex.31.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + ", paste0(vqtls[i], " * " ,vqtls[-i], collapse=" + ")))
  fit <- lm(f, pheno)
  t <- tidy(fit)
  results <- rbind(results, t[grep(":", t$term),])
}

# filter P < 0.05
results %>% filter(p.value < nrow(results) / 0.05)