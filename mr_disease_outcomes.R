library("data.table")
library("ieugwasr")
library("dplyr")
library('optparse')
library("TwoSampleMR")
set.seed(1234)

option_list = list(
  make_option(c("-i", "--id"), type="character", default=NULL, help="OpenGWAS ID", metavar="character"),
  make_option(c("-t", "--trait"), type="character", default=NULL, help="Name of trait", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, help="Output file", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# load vGWAS for biomarker risk factor
vqtl <- fread(paste0("data/", opt$trait, ".validate.txt"))

# select vQTLs with evidence of a mean effect on biomaker
vqtl <- vqtl[vqtl$Pmu < 0.05]

# select SNPs which are moderately associated using B-P
vqtl <- vqtl[vqtl$Pvar < 5e-5]

# phewas vQTLs against disease outcomes
outcomes <- phewas(vqtl$rsid, batch=c("bbj-a", "ebi-a", "finn-a", "ieu-a", "ieu-b", "ukb-a", "ukb-b", "ukb-d"))

# select outcomes which have MR evidence
mr_res <- data.frame()
for (out_id in unique(outcomes$id)){
    message(paste0("MR for ", out_id))
    # Get instruments
    exposure_dat <- extract_instruments(opt$id)

    # Get effects of instruments on outcome
    outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes=out_id)

    # Harmonise the exposure and outcome data
    dat <- harmonise_data(exposure_dat, outcome_dat)

    # Perform MR
    mr_res <- rbind(mr_res, mr(dat))
}

# select casual biomaker-outcome relationships
mr_main <- mr_res %>%
    filter(method == "Inverse variance weighted" & pval < 0.05) %>%
    pull(id.outcome)
mr_main <- as.character(mr_main)

mr_sens <- mr_res[mr_res$id.outcome %in% mr_main,] %>%
    filter(method == "Weighted median" & pval < 0.05) %>%
    pull(id.outcome)
mr_sens <- as.character(mr_sens)

# filter vQTLs with MR evidence
outcomes <- outcomes[outcomes$id %in% mr_sens,]

# write out disease outcomes
write.csv(outcomes, file=opt$out)