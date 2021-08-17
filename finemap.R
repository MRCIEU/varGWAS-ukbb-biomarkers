library("gwasglue")
library("ieugwasr")
library("susieR")
options(ieugwasr_api="http://64.227.44.193:8006/")
set.seed(13)

# lead SNP to fine map
lead <- "1:234850420"

# select natural interval around lead SNP
region <- map_variants_to_regions(chrpos=lead, pop="EUR")

# select variants and LD matrix for region
dat <- ieugwasr_to_finemapr(region$region, "ukb-d-30780_irnt", bfile = "/mnt/storage/home/ml18692/projects/jlst-cpp-vgwas/data/EUR", plink_bin = "/mnt/storage/home/ml18692/projects/jlst-cpp-vgwas/data/plink_Linux")

# perform finemapping using SuSie
fitted_rss <- susieR::susie_rss(
    dat[[1]]$z$zscore,
    dat[[1]]$ld,
    L=10,
    estimate_prior_variance=TRUE
)
summary(fitted_rss)
susieR::susie_plot(fitted_rss, y="PIP")