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

# manhattan
png(paste0("data/", opt$trait, "_manhattan.png"))
manhattan(gwas, ylim = c(0, 30), chr="CHR", bp="POS", p="P", snp="RSID", main = paste0("vGWAS Manhattan plot: ", trait_name))
dev.off()

# qq plot
png(paste0("data/", opt$trait, "_qq.png"))
qq(gwas$P, main = paste0("vGWAS Q-Q plot: ", trait_name))
dev.off()