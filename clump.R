library('data.table')
library('optparse')
library('dplyr')
library('ieugwasr')
library('broom')
source("funs.R")
set.seed(124)
options(ieugwasr_api="http://64.227.44.193:8006/")

option_list = list(
  make_option(c("-t", "--trait"), type="character", default=NULL, help="Name of trait"),
  make_option(c("-p", "--pval"), type="numeric", default=5e-8/30, help="P value threshold")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

message(paste0("trait ", opt$trait))
message(paste0("pval ", opt$p))

# load vGWAS and QC
gwas <- get_variants(opt$trait)

# select tophits
sig <- gwas[gwas$phi_p < opt$p]

# clump
stopifnot(nrow(sig)>0)
sig <- sig[,c("rsid", "phi_p")]
names(sig) <- c("rsid", "pval")
sig <- ld_clump(sig)

# subset rows
sig <- gwas[gwas$rsid %in% sig$rsid]

# flip alleles and drop effect estimate
sig <- sig %>% select(chr, pos, rsid, oa, ea) %>% rename(oa="ea", ea="oa")

# save assoc
if (opt$p == 5e-8/30){
  write.csv(sig, file=paste0("data/", opt$trait, ".clump.txt"), row.names=F)
} else {
  write.csv(sig, file=paste0("data/", opt$trait, ".clump_", opt$p,".txt"), row.names=F)
}