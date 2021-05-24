library("data.table")
library("snpStats")
set.seed(124)

# loco chr
chr <- 1

# read in snps to model
snps <- fread("/mnt/storage/private/mrcieu/research/UKBIOBANK_GWAS_Pipeline/data/model_snps_for_grm/grm6_snps.prune.in", header=F)

# read in snps in BED file
bim <- fread("/mnt/storage/private/mrcieu/research/UKBIOBANK_GWAS_Pipeline/data/bolt_bfile/grm6_european_filtered_ieu.bim")

# subset snps LOCO
bim <- bim[bim$V1 != chr]

# subset snp list to read
snps <- snps[snps$V1 %in% bim$V2]

# read in genotypes
sample <- read.plink(
    "/mnt/storage/private/mrcieu/research/UKBIOBANK_GWAS_Pipeline/data/bolt_bfile/grm6_european_filtered_ieu.bed", 
    "/mnt/storage/private/mrcieu/research/UKBIOBANK_GWAS_Pipeline/data/bolt_bfile/grm6_european_filtered_ieu.bim", 
    "/mnt/storage/private/mrcieu/research/UKBIOBANK_GWAS_Pipeline/data/bolt_bfile/grm6_european_filtered_ieu.fam",
    select.snps=snps$V1)

geno <- as(sample$genotypes, "numeric")
geno <- as.data.frame(geno)
geno$appieu <- rownames(geno)

# read in pheno
pheno <- fread("data/alanine_aminotransferase.30620.0.0.txt", select=c("appieu", "alanine_aminotransferase.30620.0.0"))

# merge
pheno <- merge(pheno, geno, "appieu")

# fit model
pheno <- na.omit(pheno)
fit <- lm(alanine_aminotransferase.30620.0.0 ~ . -appieu, data=pheno)

# get residuals
p <- predict(fit)
names(p) <- pheno$appieu

# write to file
write.table(p)