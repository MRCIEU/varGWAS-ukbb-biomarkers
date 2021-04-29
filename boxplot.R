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
  make_option(c("-c", "--chr"), type="character", default=NULL, help="Chromosome", metavar="character"),
  make_option(c("-b", "--bp"), type="character", default=NULL, help="Position", metavar="character"),
  make_option(c("-r", "--rsid"), type="character", default=NULL, help="RSID", metavar="character"),
  make_option(c("-n", "--ne"), type="character", default=NULL, help="Other allele", metavar="character"),
  make_option(c("-e", "--ea"), type="character", default=NULL, help="Effect alelle", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, help="Output file", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

message(paste0("trait ", opt$trait))

# read in extracted phenotypes
pheno <- fread(opt$p)

# GWAS
result <- model(pheno, opt$trait, as.character(opt$chr), as.numeric(opt$bp), as.character(opt$ne), as.character(opt$ea), as.character(opt$rsid))