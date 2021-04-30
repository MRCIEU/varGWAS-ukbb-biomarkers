library('data.table')
library('dplyr')
library('qqman')
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

# load vGWAS for biomarker risk factor
gwas <- data.frame()
for (chr in seq(16,22)){
    if (chr < 10){
        file <- paste0("data/", opt$trait, ".vgwas.chr0", chr, ".txt")
    } else {
        file <- paste0("data/", opt$trait, ".vgwas.chr", chr, ".txt")
    }
    gwas <- rbind(gwas, fread(file))
}

# drop HLA region
gwas <- gwas[!(gwas$chr == 6 & gwas$pos >= 29691116 & gwas$pos <= 33054976),]

# drop failed rows
gwas <- gwas %>%
    filter(n != -1)

# manhattan
png(paste0("data/", opt$trait, "_mean_manhattan.png"))
manhattan(gwas, ylim = c(0, 25), chr="chr", bp="pos", p="p", snp="rsid", main = trait_name)
dev.off()

png(paste0("data/", opt$trait, "_phi_manhattan.png"))
manhattan(gwas, ylim = c(0, 25), chr="chr", bp="pos", p="phi_p", snp="rsid", main = trait_name)
dev.off()

# qq plot
png(paste0("data/", opt$trait, "_mean_qq.png"))
qq(gwas$p, main = trait_name)
dev.off()

png(paste0("data/", opt$trait, "_phi_qq.png"))
qq(gwas$phi_p, main = trait_name)
dev.off()