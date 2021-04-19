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

# select vQTLs with evidence of a mean effect on outcome
vqtl <- vqtl[vqtl$Pmu < 0.01]

# select SNPs which are strongly associated with B-P
vqtl <- vqtl[vqtl$Pvar < 5e-8]

# phewas vQTL
exposures <- phewas(vqtl$rsid)

# drop traits with high genetic correlation with outcome
# TODO

# select exposures associated with multiple vQTLs
# TODO optimise threshold
exposures.ms <- exposures[exposures$id %in% exposures$id[table(exposures$id) > 3],]

# colocalise vQTL with exposure GWAS to check for shared causal variant
# TODO

# output candidate modifiers for vQTLs
# TODO

# MVMR & observational interaction effect of exposure & modifier on outcome
# TODO