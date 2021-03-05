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

mod <- function(pheno, out, chr, pos, oa, ea) {
  tryCatch(
    {
      dosage <- extract_variant_from_bgen(chr, pos, oa, ea)
      pheno <- merge(pheno, dosage, "appieu")
      pheno <- na.omit(pheno)
      s <- paste0("chr", chr, "_", pos, "_", oa, "_", ea)
      f <- paste0(out, " ~ ", s, " + sex.31.0.0 + age_at_recruitment.21022.0.0 +", paste0("PC", seq(1, 10), collapse="+"))
      fit1 <- lm(as.formula(f), data=pheno)
      pheno$d <- resid(fit1)^2
      f <- paste0("d ~ ", s, " + sex.31.0.0 + age_at_recruitment.21022.0.0 +", paste0("PC", seq(1, 10), collapse="+"))
      fit2 <- lm(as.formula(f), data=pheno)
      return(tidy(fit2)[2,])
  },
  error=function(cond) {
      return(NA)
    }
  )
}

# read in extracted phenotypes
pheno <- fread(opt$p)

# read in snp list
snps <- fread(opt$s)

# GWAS
results <- apply(snps, 1, function(snp) {
  mod(pheno, opt$t, as.character(snp[['chromosome']]), as.numeric(snp[['position']]), snp[['first_allele']], snp[['alternative_alleles']])
})
results <- rbindlist(results)

# save assoc
write.csv(results, file=opt$o, row.names=F)