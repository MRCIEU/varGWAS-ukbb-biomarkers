library('data.table')
library('optparse')
library('dplyr')
library('ieugwasr')
library('quantreg')
library('broom')
library("robustbase")
source("funs.R")
source("model.R")
set.seed(124)

option_list = list(
  make_option(c("-p", "--pheno_file"), type="character", default=NULL, help="UKBB pheno CSV", metavar="character"),
  make_option(c("-t", "--trait"), type="character", default=NULL, help="Name of trait", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, help="Output file", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

message(paste0("trait ", opt$trait))

# read in extracted phenotypes
pheno <- fread(opt$p)

# load vGWAS and QC
gwas <- get_variants(opt$trait)

# select tophits
sig <- gwas[gwas$P < (5e-8/30)]

# clump
stopifnot(nrow(sig)>0)
sig <- sig[,c("RSID", "P")]
names(sig) <- c("rsid", "pval")
sig <- ld_clump(sig)

# subset rows
sig <- gwas[gwas$RSID %in% sig$rsid]

# drop multialllelic variants
ma <- as.data.frame(table(sig$RSID)) %>% filter(Freq == 2) %>% mutate_if(is.factor, as.character) %>% pull(Var1)
sig <- sig[!(sig$RSID %in% ma)]

# GWAS
message(paste0("Found ", nrow(sig), " variants"))
results <- apply(sig, 1, function(snp) {
  model(pheno, opt$trait, as.character(snp[['CHR']]), as.numeric(snp[['POS']]), as.character(snp[['OA']]), as.character(snp[['EA']]), as.character(snp[['RSID']]))
})
results <- rbindlist(results[!is.na(results)], fill=T)

# save assoc
write.csv(results, file=opt$o, row.names=F)