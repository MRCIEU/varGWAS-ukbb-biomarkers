library("ieugwasr")
source("funs.R")
set.seed(1234)
options(ieugwasr_api="http://64.227.44.193:8006/")

# load BMI vGWAS
vgwas <- get_variants("body_mass_index.21001.0.0")

# subset GLP1R
GLP1R <- vgwas %>% dplyr::filter(position >= 39016574-500000 & position <= 39055519+500000)

# clump
p <- GLP1R %>% dplyr::filter(position >= 39016574-500000 & position <= 39055519+500000) %>% dplyr::select(rsid, phi_p) %>% dplyr::rename(p="phi_p")
vqtl <- ieugwasr::ld_clump(p)

# extract QTL tophits
qtl <- tophits(c("prot-a-1219", "finn-a-GLP1ANA", "eqtl-a-ENSG00000112164"))

# test for vQTL effect
