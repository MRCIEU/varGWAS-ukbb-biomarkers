library("ieugwasr")
library("dplyr")
library("metafor")
source("funs.R")
set.seed(1234)

get_rct_est <- function(placebo_mean, placebo_n, placebo_sd, treatment_mean, treatment_n, treatment_sd){
    # Taken from https://github.com/harrietlmills/DetectingDifferencesInVariance/blob/master/AnalyseIndividualTrials.R
    rdat_logVR <- escalc(measure = "VR", m1i = treatment_mean, n1i = treatment_n, sd1i = treatment_sd,  m2i = placebo_mean, n2i = placebo_n, sd2i = placebo_sd)

    # calculate confidence intervals, etc
    srdat_logVR <- summary(rdat_logVR, digits = 2)

    # calculate test statistic and pvalue
    logVR_pvalue <-  2*(1-pnorm(abs(srdat_logVR$zi)))

    est <- rdat_logVR$yi 
    se <-  srdat_logVR$sei 
    test <- srdat_logVR$zi
    p <- logVR_pvalue

    return(est, se, test, p)
}

# load instruments for lipids at drug target loci & cross ref with vGWAS

# LDL-c

# get instruments
ldl <- tophits("ieu-b-110")

# load vGWAS for LDL-c
vgwas <- get_variants("ldl_direct.30780.0.0")

# filter drug target loci
hmgcr <- ldl %>% filter(chr == "5" & position > 74632154-500000 & position < 74657929+500000)
hmgcr <- vgwas %>% filter(rsid %in% hmgcr$rsid)
### PCSK9 effect on LDL-c ###
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4845239/
# RCT - PCSK9 inhibitors lower LDL-c mean and variance (p=0.005; n=29)
# MR - rs191448950:A at PCSK9 lowers LDL-c mean and variance (p=3.15 x 10-8)
pcsk9 <- ldl %>% filter(chr == "1" & position > 55505221-500000 & position < 55530525+500000)
pcsk9 <- vgwas %>% filter(rsid %in% pcsk9$rsid)
npc1l1 <- ldl %>% filter(chr == "7" & position > 44552134-500000 & position < 44580914+500000)
npc1l1 <- vgwas %>% filter(rsid %in% npc1l1$rsid)
cetp <- ldl %>% filter(chr == "16" & position > 56995762-500000 & position < 57017757+500000)
cetp <- vgwas %>% filter(rsid %in% cetp$rsid)

# HDL-c

# get instruments
hdl <- tophits("ieu-b-109")

# load vGWAS for HDL-c
vgwas <- get_variants("hdl_cholesterol.30760.0.0")

# filter drug target loci
cetp <- hdl %>% filter(chr == "16" & position > 56995762-500000 & position < 57017757+500000)
cetp <- vgwas %>% filter(rsid %in% cetp$rsid)

### vGWAS shows reduction in variance at drug target loci, does RCT show same?

# RCT evidence
# TODO
get_rct_est()