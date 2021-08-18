library("gwasglue")
library("ieugwasr")
library("data.table")
library('optparse')
library("stringr")
library("susieR")
options(ieugwasr_api="http://64.227.44.193:8006/")
set.seed(13)

option_list = list(
  make_option(c("-t", "--trait"), type="character", default=NULL, help="Name of trait")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# load clumped vQTLs
d <- fread("data/vqtls.txt")
d <- d %>% dplyr::filter(trait != "body_mass_index.21001.0.0")
d$key <- paste0("chr", d$chr, "_", d$pos, "_", d$oa, "_", d$ea)
d$id <- paste0("ukb-d-", str_split(d$trait, "\\.", simplify=T)[,2], "_irnt")
d$chr_pos <- paste0(d$chr, ":", d$pos)

for (i in 1:nrow(d)){
    # select natural interval around lead SNP
    region <- map_variants_to_regions(chrpos=d$chr_pos[i], pop="EUR")

    # select variants and LD matrix for region
    dat <- ieugwasr_to_finemapr(
        region$region,
        d$id[i],
        bfile = "/mnt/storage/home/ml18692/projects/jlst-cpp-vgwas/data/EUR",
        plink_bin = "/mnt/storage/home/ml18692/projects/jlst-cpp-vgwas/data/plink_Linux"
    )

    # perform finemapping using SuSie
    fitted_rss <- susieR::susie_rss(
        dat[[1]]$z$zscore,
        dat[[1]]$ld,
        L=10,
        estimate_prior_variance=TRUE
    )

    # collect fine mapped snps
    snps <- character()
    for (i in fitted_rss$sets$cs) { 
        for (j in i) {
            snps <- c(snps, dat[[1]]$z$snp[j])
        }
    }

    return(associations(snps, id))
}