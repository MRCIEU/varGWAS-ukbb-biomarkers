library("data.table")
library("ieugwasr")
set.seed(1234)

trait <- "alanine_aminotransferase.30620"

vqtl <- data.frame()
for (chr in seq(1,22)){
    if (chr < 10){
        file <- paste0("data/", trait, ".validate.chr0", chr, ".txt")
    } else {
        file <- paste0("data/", trait, ".validate.chr", chr, ".txt")
    }
    vqtl <- rbind(vqtl, fread(file))
}

# select SNPs which are strongly associated with B-P
vqtl <- vqtl[vqtl$P < 5e-8]

# phewas vQTL
exposures <- phewas(vqtl$rsid)

# select exposures with multiple SNPs
exposures.ms <- exposures[exposures$id %in% exposures$id[table(exposures$id) > 3],]