library('data.table')
library('optparse')
library('dplyr')
library('ieugwasr')
library('broom')
source("funs.R")
set.seed(124)

option_list = list(
  make_option(c("-t", "--trait"), type="character", default=NULL, help="Name of trait", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

message(paste0("trait ", opt$trait))

# load vGWAS and QC
gwas <- get_variants(opt$trait)

# select tophits
sig <- gwas[gwas$phi_p < (5e-8/30)]
sig <- sig[sig$p < (0.1/nrow(sig))]

# clump
stopifnot(nrow(sig)>0)
sig <- sig[,c("rsid", "phi_p")]
names(sig) <- c("rsid", "pval")
sig <- ld_clump(sig)

# subset rows
sig <- gwas[gwas$rsid %in% sig$rsid]
sig$key <- NULL

# save assoc
write.csv(sig, file=paste0("data/", opt$trait, ".clump.txt"), row.names=F)