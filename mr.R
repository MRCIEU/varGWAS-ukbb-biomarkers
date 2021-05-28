library("data.table")
library("TwoSampleMR")
source("funs.R")
set.seed(1234)

option_list = list(
  make_option(c("-c", "--trait"), type="character", default=NULL, help="Name of trait", metavar="character"),
  make_option(c("-o", "--opengwas"), type="character", default=NULL, help="biomarker ID in opengwas", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# load vGWAS and QC
data <- get_variants(opt$trait)

# read in coloc data
coloc <- fread(paste0("data/", opt$trait, ".coloc.txt"))

# retain high prob of shared casual variant
coloc <- coloc[coloc$PP.H4.abf > .8]

results <- data.frame()
for (id in unique(coloc$gene)){
    message(paste0("Working on id: ", id))

    # Test for effect of gene product on biomarker conc

    # Get instruments
    exposure_dat <- extract_instruments(id)

    # Get effects of instruments on outcome
    outcome_dat <- data[data$rsid %in% exposure_dat$SNP]
    outcome_dat <- format_data(
        outcome_dat,
        type="outcome",
        snp_col="rsid",
        beta_col="beta", 
        se_col="se",
        eaf_col="eaf",
        effect_allele_col="ea",
        other_allele_col="oa", 
        pval_col="p",
        samplesize_col="n"
    )

    # Harmonise the exposure and outcome data
    dat <- harmonise_data(exposure_dat, outcome_dat)

    # Perform MR
    res <- mr(dat)
    res <- cbind(res, directionality_test(dat))
    results <- rbind(results, res)

    # Test for effect of biomarker on gene product

    # Get instruments
    exposure_dat <- extract_instruments(opt$o)

    # Get effects of instruments on outcome
    outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes=id)

    # Harmonise the exposure and outcome data
    dat <- harmonise_data(exposure_dat, outcome_dat)

    # Perform MR
    res <- mr(dat)
    res <- cbind(res, directionality_test(dat))
    results <- rbind(results, res)
}