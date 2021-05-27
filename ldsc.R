library('data.table')
library('dplyr')
library('qqman')
library('stringr')
library('optparse')
source("funs.R")
set.seed(1234)

option_list = list(
  make_option(c("-t", "--trait"), type="character", default=NULL, help="Name of trait", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

message(paste0("trait ", opt$trait))

# load vGWAS and QC
data <- get_variants(opt$trait)

# stringent filter on INFO
data <- data[which(data$info > .9),]

# select fields for LDSC
write.table(data[,c("rsid", "ea", "oa", "n", "p", "beta", "eaf", "info")], sep="\t", quote=F, row.names=F, file=paste0("data/", opt$trait, ".ldsc"))