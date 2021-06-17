library('data.table')
library('dplyr')
library('stringr')
library('optparse')
library('ggplot2')
library('car')
source("funs.R")
set.seed(1234)

hists <- list()
qqs <- list()

for (i in 1:length(biomarkers)){
    trait_name <- get_trait_name(biomarkers[i])
    labs[i] <- trait_name
    message(paste0("trait ", biomarkers[i]))
    message(paste0("trait name ", trait_name))

    # load phenotype
    data <- fread(paste0("data/", opt$trait, ".txt"))

    # plots
    hists[[i]] <- hist(d[[opt$trait]], xlab=trait_name, main=trait_name)
    qqs[[i]] <- qqPlot(d[[opt$trait]], main = trait_name, ylab=trait_name)
}

option_list = list(
  make_option(c("-t", "--trait"), type="character", default=NULL, help="Name of trait", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

trait_name <- get_trait_name(opt$trait)
message(paste0("trait ", opt$trait))
message(paste0("trait name ", trait_name))

# load phenotype
d <- fread(paste0("data/", opt$trait, ".txt"))

# histogram
png(paste0("data/", opt$trait, "_pheno_hist.png"))
hist(d[[opt$trait]], xlab=trait_name, main=trait_name)
dev.off()

# qq plot
png(paste0("data/", opt$trait, "_pheno_qq.png"))
qqPlot(d[[opt$trait]], main = trait_name, ylab=trait_name)
dev.off()