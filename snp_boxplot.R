library('data.table')
library('optparse')
library('dplyr')
library('ieugwasr')
library('quantreg')
library('broom')
library("robustbase")
library("scales")
library("ggplot2")
source("funs.R")
set.seed(124)

option_list = list(
  make_option(c("-p", "--pheno_file"), type="character", default=NULL, help="UKBB pheno CSV", metavar="character"),
  make_option(c("-t", "--trait"), type="character", default=NULL, help="Name of trait", metavar="character"),
  make_option(c("-c", "--chr"), type="character", default=NULL, help="Chromosome", metavar="character"),
  make_option(c("-b", "--bp"), type="character", default=NULL, help="Position", metavar="character"),
  make_option(c("-n", "--ne"), type="character", default=NULL, help="Other allele", metavar="character"),
  make_option(c("-e", "--ea"), type="character", default=NULL, help="Effect alelle", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, help="Output file", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

opt <- data.frame(trait="alanine_aminotransferase.30620", chr="4", bp="88212722", ne="G", ea="A", p="data/alanine_aminotransferase.30620.txt", stringsAsFactors=F)

message(paste0("trait ", opt$trait))

# read in extracted phenotypes
data <- fread(opt$p)

# get SNP
dosage <- extract_variant_from_bgen(as.character(opt$chr), as.numeric(opt$bp), as.character(opt$ne), as.character(opt$ea))

# prepare data
data <- merge(data, dosage, "appieu")
data <- na.omit(data)
s <- paste0("chr", as.character(opt$chr), "_",  as.numeric(opt$bp), "_", as.character(opt$ne), "_", as.character(opt$ea))
s <- gsub(" ", "", s, fixed = TRUE)
data[[s]] <- as.factor(round(data[[s]]))

# produce plot
p <- ggplot(data, aes_string(x=s, y=opt$trait)) + geom_boxplot() + theme_classic() + scale_y_continuous(trans = log2_trans(), breaks = trans_breaks("log2", function(x) 2^x), labels = trans_format("log2", math_format(2^.x)))

# mean and variance by genotype group
data %>%
    group_by(chr4_88212722_G_A) %>%
    summarise(mean = mean(alanine_aminotransferase.30620), var = var(alanine_aminotransferase.30620))