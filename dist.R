library('data.table')
library('dplyr')
library('car')
library('stringr')
library('optparse')
set.seed(1234)

option_list = list(
  make_option(c("-t", "--trait"), type="character", default=NULL, help="Name of trait", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

trait_name <- opt$trait
trait_name <- str_split(trait_name, "\\.", simplify = TRUE)[,1]
trait_name <- gsub("_", " ", trait_name)
trait_name <- str_to_title(trait_name)
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