library('data.table')
library('optparse')
library('dplyr')
library('ieugwasr')
library('quantreg')
library('broom')
source("funs.R")
set.seed(124)

option_list = list(
  make_option(c("-p", "--pheno_file"), type="character", default=NULL, help="UKBB pheno CSV", metavar="character"),
  make_option(c("-t", "--trait"), type="character", default=NULL, help="Name of trait", metavar="character"),
  make_option(c("-g", "--gwas"), type="character", default=NULL, help="Path to GWAS summary stats", metavar="character"),
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
      s2 <- paste0(s, "2")
      pheno[[s2]] <- pheno[[s]]^2
      f <- paste0(out, " ~ ", s, " + sex.31.0.0 + age_at_recruitment.21022.0.0 +", paste0("PC", seq(1, 10), collapse="+"))
      fit1 <- rq(as.formula(f), tau=0.5, data=pheno)
      pheno$d <- resid(fit1)^2
      f <- paste0("d ~ ", s, " + ", s2)
      fit2 <- lm(as.formula(f), data=pheno)
      f <- paste0("d ~ 1")
      fit0 <- lm(as.formula(f), data=pheno)
      ftest <- tidy(anova(fit0, fit2))
      fit2 <- tidy(fit2)
      return(list(ftest, fit2))
    },
    error=function(cond) {
      return(NA)
    }
  )
}

# read in extracted phenotypes
pheno <- fread(opt$p)

# read in GWAS
gwas <- fread(opt$g)

# drop failed rows
gwas <- gwas %>%
    filter(SE_x != -1)

# select tophits
sig <- gwas[gwas$P < 5e-8 & gwas$EAF >= 0.03 & gwas$EAF <= 0.98]
sig <- sig[,c("RSID", "P")]
names(sig) <- c("rsid", "pval")
sig <- ld_clump(sig)

# subset rows
sig <- gwas[gwas$RSID %in% sig$rsid]

# GWAS
results <- apply(sig, 1, function(snp) {
  mod(pheno, opt$t, as.character(snp[['CHR']]), as.numeric(snp[['POS']]), as.character(snp[['OA']]), as.character(snp[['EA']]))
})
results <- rbindlist(results[!is.na(results)], fill=T)

# save assoc
write.csv(results, file=opt$o, row.names=F)