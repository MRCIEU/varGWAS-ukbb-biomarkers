library("data.table")
library("ieugwasr")
library("dplyr")
library("TwoSampleMR")
set.seed(1234)

trait <- "calcium.30680"
trait_id <- "ukb-d-30680_irnt"

vqtl <- data.frame()
for (chr in seq(1,22)){
    if (chr < 10){
        file <- paste0("data/", trait, ".validate.chr0", chr, ".txt")
    } else {
        file <- paste0("data/", trait, ".validate.chr", chr, ".txt")
    }
    if (file.exists(file)){
        vqtl <- rbind(vqtl, fread(file))
    } else {
        warning(paste0("File ", file, " does not exist"))
    }
}

# select vQTLs with evidence of a mean effect on outcome
vqtl <- vqtl[vqtl$Pmu < 0.05]

# select SNPs which are strongly associated with B-P
vqtl <- vqtl[vqtl$Pvar < 5e-5]

# phewas vQTL
exposures <- phewas(vqtl$rsid)

# select exposures associated with multiple vQTLs
# TODO optimise threshold
counts <- as.data.frame(table(exposures$id), stringsAsFactors=F)
ids <- counts[counts$Freq > 1,]$Var1
exposures.ms <- exposures[exposures$id %in% ids,]

# MR of exposures/modifiers on outcome
mr_res <- data.frame()
for (exp_id in unique(exposures.ms$id)){
    message(paste0("MR for ", exp_id))
    # Get instruments
    exposure_dat <- extract_instruments(exp_id)

    # Get effects of instruments on outcome
    outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes=trait_id)

    # Harmonise the exposure and outcome data
    dat <- harmonise_data(exposure_dat, outcome_dat)

    # Perform MR
    mr_res <- rbind(mr_res, mr(dat), stringsAsFactors=F)
}

# select exposures with a causal effect on outcome
mr_exp_effect <- mr_res %>%
    filter(method == "Inverse variance weighted" & pval < 0.05) %>%
    pull(id.exposure)
mr_exp_effect <- as.character(mr_exp_effect)

# MR sensitivity analysis
mr_exp_effect <- mr_res[mr_res$id.exposure %in% mr_exp_effect,] %>%
    filter(method == "Weighted median" & pval < 0.05) %>%
    pull(id.exposure)
mr_exp_effect <- as.character(mr_exp_effect)

# filter vqtls where there is MR evidence of an exposure/modifier effect on outcome
exposures.mr <- exposures.ms[exposures.ms$id %in% mr_exp_effect,]