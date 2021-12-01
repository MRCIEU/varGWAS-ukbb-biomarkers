library("ieugwasr")
library("ggplot2")
library("viridis")
library("RColorBrewer")
source("funs.R")
set.seed(123)

# region
chr <- "4"
start <- 10402838 - 500000
end <- 10402838 + 500000

# load vGWAS
vgwas <- get_variants("urate.30880.0.0")
vqtl <- vgwas %>%
    dplyr::filter(chr == !!chr) %>% 
    dplyr::filter(pos >= !!start & pos <= !!end)

# load LD
ld_dat <- ieugwasr::ld_matrix (
    vqtl$rsid,
    with_alleles = FALSE,
    pop = "EUR",
    bfile = "/user/home/ml18692/varGWAS-ukbb-biomarkers/data/EUR",
    plink_bin = "/user/home/ml18692/varGWAS-ukbb-biomarkers/data/plink_Linux"
)

# add LD for top vQTL (rs4530622)
vqtl <- vqtl %>% dplyr::filter(rsid %in% row.names(ld_dat))
vqtl$ld <- ld_dat[row.names(ld_dat) == "rs4530622"]^2
vqtl <- vqtl %>% filter(ld > 0.05)
vqtl$lp <- -log10(vqtl$p)
vqtl$lvp <- -log10(vqtl$phi_p)
vqtl <- vqtl %>% mutate(ld1 = 
    case_when(
        ld > .8 ~ 4,
        ld > .6 ~ 3,
        ld > .4 ~ 2,
        ld > .2 ~ 1,
        TRUE ~ 0
    )
)
vqtl$ld1 <- as.factor(vqtl$ld1)
ggplot(vqtl, aes(x=lp, y=lvp, color=ld1)) +
    geom_point() +
    scale_color_manual(values=c("darkblue", "lightblue", "green", "orange", "red")) +
    geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "grey") +
    geom_vline(xintercept = -log10(5e-8), linetype = "dashed", color = "grey") +
    theme_bw()
