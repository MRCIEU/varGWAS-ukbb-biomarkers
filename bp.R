library('optparse')
library('data.table')
library('broom')
source("funs.R")
set.seed(1234)

# Run B-P model using R to validate cpp application

option_list = list(
  make_option(c("-p", "--pheno_file"), type="character", default=NULL, help="UKBB pheno CSV", metavar="character"),
  make_option(c("-t", "--trait"), type="character", default=NULL, help="Name of trait", metavar="character"),
  make_option(c("-s", "--snp_file"), type="character", default=NULL, help="SNP file", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, help="Output file", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

opt <- data.frame(p="data/ukb_bmi.txt", t="body_mass_index.21001.0.0", s="data/snps.txt", o="data/ukb_bmi.vgwas.r_subsample.txt", stringsAsFactors=F)

mod <- function(pheno, out, chr, pos, oa, ea) {
    dosage <- extract_variant_from_bgen(chr, pos, oa, ea)
    pheno <- merge(pheno, dosage, "appieu")
    fit1 <- lm(pheno[[out]] ~ )
    return(tidy(fit)[2,])
}

# read in extracted phenotypes
pheno <- fread(opt$p)

# read in snp list
snps <- fread(opt$s)
snp <- snps[1]

# GWAS
results <- apply(snps, 1, function(snp) {
  mod(pheno, opt$t, snp[['chr']], as.numeric(snp[['pos']]), snp[['other_allele']], snp[['effect_allele']])
})
results <- rbindlist(results)

# save assoc
write.csv(results, file=opt$o, row.names=F)