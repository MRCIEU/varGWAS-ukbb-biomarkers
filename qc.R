library('data.table')
library('dplyr')
library('qqman')
library('optparse')
set.seed(1234)

option_list = list(
  make_option(c("-t", "--trait"), type="character", default=NULL, help="Name of trait", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

message(paste0("trait ", opt$trait))

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
gwas$Ptemp <- gwas$P
gwas$Ptemp[gwas$Ptemp == 0] <- .Machine$double.xmin
png(paste0("data/", opt$trait, "_manhattan.png"))
manhattan(gwas, chr="CHR", bp="POS", p="Ptemp", snp="RSID", main = paste0("Manhattan plot of variance GWAS p-values: ", opt$trait))
dev.off()

# qq plot
png(paste0("data/", opt$trait, "_qq.png"))
qq(gwas$P, main = paste0("Q-Q plot of variance GWAS p-values: ", opt$trait))
dev.off()