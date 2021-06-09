library('optparse')
library('data.table')
library('broom')
library('jlst')
source("funs.R")
set.seed(1234)

# Run B-P model using R to validate cpp application

option_list = list(
  make_option(c("-t", "--trait"), type="character", default=NULL, help="Name of trait", metavar="character"),
  make_option(c("-m", "--model"), type="character", default=NULL, help="Model", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

mod <- function(pheno, out, chr, pos, oa, ea, model) {
  tryCatch(
    {
      dosage <- extract_variant_from_bgen(chr, pos, oa, ea)
      pheno <- merge(pheno, dosage, "appieu")
      pheno <- na.omit(pheno)
      s <- paste0("chr", chr, "_", pos, "_", oa, "_", ea)

      if (model == "BP"){
        t <- vartest(pheno[[out]], pheno[[s]], covar=pheno[,c("age_at_recruitment.21022.0.0", "sex.31.0.0", paste0("PC", seq(1, 10)))], covar.var=T, type=1, x.sq=T)
        return(data.frame(p=t$test$P, snp=s, model))
      } else if (model == "JLSSC"){
        t <- jlssc(pheno[[out]], pheno[[s]], covar=pheno[,c("age_at_recruitment.21022.0.0", "sex.31.0.0", paste0("PC", seq(1, 10)))], type=1, x.sq=T)
        return(data.frame(p=t$P, snp=s, model))
      } else if (model == "MOM"){
        t <- jlssc(pheno[[out]], pheno[[s]], covar=pheno[,c("age_at_recruitment.21022.0.0", "sex.31.0.0", paste0("PC", seq(1, 10)))], type=3, x.sq=T)
        return(data.frame(p=t$P, snp=s, model))
      }

    },
    error=function(cond) {
      return(NA)
    }
  )
}

message(paste0("vgwas of ", opt$t))

# read in extracted phenotypes
pheno <- fread(paste0("data/", opt$t, ".txt"))

# read in snp list
snps <- fread(paste0("data/", opt$t, ".30k_snps.txt"))

# GWAS
results <- apply(snps, 1, function(snp) {
  mod(pheno, opt$t, as.character(snp[['chromosome']]), as.numeric(snp[['position']]), as.character(snp[['first_allele']]), as.character(snp[['alternative_alleles']]), opt$m)
})
results <- rbindlist(results[!is.na(results)])

# save assoc
write.table(results, file=paste0("data/", opt$t, ".30k_snps.", opt$m,".txt"), row.names=F,quote=F, sep="\t")