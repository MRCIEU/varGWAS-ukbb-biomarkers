library("coloc")
library("ieugwasr")
source("funs.R")
set.seed(123)

# region
chr <- "4"
start <- 10402838 - 1250000
end <- 10402838 + 500000

# load vGWAS
vgwas <- get_variants("urate.30880.0.0")

# extract 1Mb interval around vQTL
qtl <- vgwas %>%
    dplyr::filter(chr == !!chr) %>% 
    dplyr::filter(pos >= !!start & pos <= !!end) %>%
    dplyr::filter(grepl("^rs", rsid)) %>%
    dplyr::mutate(varbeta = se^2) %>% 
    dplyr::select(beta, varbeta, rsid, oa, ea, pos, n, minor_allele_frequency, p) %>%
    dplyr::rename(snp="rsid", position="pos", N="n", MAF="minor_allele_frequency", pvalues="p")

vqtl <- vgwas %>%
    dplyr::filter(chr == !!chr) %>% 
    dplyr::filter(pos >= !!start & pos <= !!end) %>%
    dplyr::filter(grepl("^rs", rsid)) %>%
    dplyr::select(rsid, oa, ea, pos, n, minor_allele_frequency, phi_p) %>%
    dplyr::rename(snp="rsid", position="pos", N="n", MAF="minor_allele_frequency", pvalues="phi_p")

# select variants and LD matrix for region
stopifnot(all(qtl$snp == vqtl$snp))
ld_dat <- ieugwasr::ld_matrix (
    qtl$snp,
    with_alleles = TRUE,
    pop = "EUR",
    bfile = "/user/home/ml18692/varGWAS-ukbb-biomarkers/data/EUR",
    plink_bin = "/user/home/ml18692/varGWAS-ukbb-biomarkers/data/plink_Linux"
)

qtl$snp <- paste0(qtl$snp, "_", qtl$oa, "_", qtl$ea)
vqtl$snp <- paste0(vqtl$snp, "_", vqtl$oa, "_", vqtl$ea)

qtl <- qtl %>% dplyr::filter(snp %in% row.names(ld_dat))
vqtl <- vqtl %>% dplyr::filter(snp %in% row.names(ld_dat))

# format as coloc data
qtl <- as.list(qtl)
qtl$type <- "quant"
qtl$LD <- ld_dat
qtl$sdY <- 1
check_dataset(qtl,warn.minp=1e-10)
#plot_dataset(qtl)

vqtl <- as.list(vqtl)
vqtl$type <- "quant"
vqtl$LD <- ld_dat
vqtl$sdY <- 1
check_dataset(vqtl,warn.minp=1e-10)
#plot(-log10(vqtl$pvalues))

# SuSiE
S_vqtl=runsusie(vqtl)
S_qtl=runsusie(qtl)