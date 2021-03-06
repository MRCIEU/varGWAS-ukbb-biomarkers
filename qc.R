library('data.table')
library('dplyr')
library('qqman')

# read in GWAS data
gwas <- data.frame()
for (chr in seq(1,22)){
    if (chr < 10){
        path <- paste0("data/ukb_bmi.vgwas.chr0", chr, ".txt")
    } else {
        path <- paste0("data/ukb_bmi.vgwas.chr", chr, ".txt")
    }
    message(paste0("Reading GWAS: ", path))
    gwas <- rbind(gwas, fread(path))
}

# drop failed rows
gwas <- gwas %>%
    filter(SE != -1 & P != -1)

# manhattan
png("manhattan.png")
manhattan(gwas, chr="CHR", bp="POS", p="P", snp="RSID", main = "Manhattan plot of variance GWAS p-values")
dev.off()

# qq plot
png("qq.png")
qq(gwas$P, main = "Q-Q plot of variance GWAS p-values")
dev.off()