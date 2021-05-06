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

trait_name <- opt$trait
trait_name <- str_split(trait_name, "\\.", simplify = TRUE)[,1]
trait_name <- gsub("_", " ", trait_name)
trait_name <- str_to_title(trait_name)
message(paste0("trait ", opt$trait))
message(paste0("trait name ", trait_name))

# load vGWAS and QC
data <- get_variants(opt$trait)

# manhattan
png(paste0("data/", opt$trait, "_phi_manhattan.png"))
manhattan(data, ylim = c(0, 25), chr="chr", bp="pos", p="phi_p", snp="rsid", main = trait_name)
dev.off()

# qq plot
png(paste0("data/", opt$trait, "_phi_qq.png"))
qq(data$phi_p, main = trait_name)
dev.off()