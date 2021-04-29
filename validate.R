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

# load vGWAS for biomarker risk factor
gwas <- data.frame()
for (chr in seq(1,22)){
    if (chr < 10){
        file <- paste0("data/", opt$trait, ".vgwas.chr0", chr, ".txt")
    } else {
        file <- paste0("data/", opt$trait, ".vgwas.chr", chr, ".txt")
    }
    gwas <- rbind(gwas, fread(file))
}

# drop HLA region
gwas <- gwas[!(gwas$CHR == 6 & gwas$POS >= 29691116 & gwas$POS <= 33054976),]

# drop failed rows
gwas <- gwas %>%
    filter(SE_x != -1)

# select tophits
sig <- gwas[gwas$P < (5e-8/30) & gwas$EAF >= 0.05 & gwas$EAF <= 0.95]

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