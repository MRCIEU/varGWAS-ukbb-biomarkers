library('data.table')
library('optparse')
library('dplyr')
library('ieugwasr')
library('jlst')
library('robustbase')
library('broom')
source("funs.R")
set.seed(124)

option_list = list(
  make_option(c("-t", "--trait"), type="character", default=NULL, help="Name of trait", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

model <- function(data, out, chr, pos, oa, ea) {
  # prepare data
  dosage <- extract_variant_from_bgen(chr, pos, oa, ea)
  data <- merge(data, dosage, "appieu")
  data <- na.omit(data)
  s <- paste0("chr", chr, "_", pos, "_", oa, "_", ea)
  s <- gsub(" ", "", s, fixed = TRUE)

  # test for mean effect using robust model
  f <- as.formula(paste0(out, " ~ ", s, " + sex.31.0.0 + age_at_recruitment.21022.0.0 +", paste0("PC", seq(1, 10), collapse="+")))
  fit <- tidy(lmrob(f, data=data))[2,]
  names(fit) <- c("key", "beta.robust", "se.robust", "t.robust", "p.robust")

  # test for mean-variance effect using B-P method of moments
  j <- jlssc(data[[out]], data[[s]], covar=data[,c("sex.31.0.0", "age_at_recruitment.21022.0.0", paste0("PC", seq(1, 10)))], type=3, x.sq = T)
  names(j)<-c("Q.jlssc", "DF.jlssc", "P.jlssc")

  return(cbind(j, fit))
}

message(paste0("trait ", opt$trait))

# load vGWAS and QC
gwas <- get_variants(opt$trait)

# select tophits
sig <- gwas[gwas$phi_p < (5e-8/30)]

# clump
stopifnot(nrow(sig)>0)
sig <- sig[,c("rsid", "phi_p")]
names(sig) <- c("rsid", "pval")
sig <- ld_clump(sig)

# subset rows
sig <- gwas[gwas$rsid %in% sig$rsid]

# test for mean effect on trait using heteroscedastiy robust model & mean-variance effect using method of moments

# read in extracted phenotypes
pheno <- fread(paste0("data/", opt$trait, ".txt"))

results <- apply(sig, 1, function(snp) {
  model(pheno, opt$trait, as.character(snp[['chr']]), as.numeric(snp[['pos']]), as.character(snp[['oa']]), as.character(snp[['ea']]))
})
results <- rbindlist(results[!is.na(results)], fill=T)

# merge results with clumped data
sig$key <- paste0("chr", sig$key)
sig <- merge(sig, results, "key")

# filter out vQTLs that do not have a mean effect
sig <- sig %>% filter(P.jlssc < 5e-5 & p.robust < 5e-5)

# save assoc
sig$key <- NULL
write.csv(sig, file=paste0("data/", opt$trait, ".clump.txt"), row.names=F)