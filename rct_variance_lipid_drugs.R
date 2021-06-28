library("metafor")
library("readxl")
library('meta')
library("dplyr")
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

    return(data.frame(est, se, test, p))
}

# RCT 
dat <- read_excel("RCT.xlsx")
dat$mu_placebo <- as.numeric(dat$mu_placebo)
dat$mu_treatment <- as.numeric(dat$mu_treatment)
dat$sigma_placebo <- as.numeric(dat$sigma_placebo)
dat$sigma_treatment <- as.numeric(dat$sigma_treatment)
dat <- na.omit(dat)

# remove study which has no placebo
dat <- dat %>% filter(Study != "ILLUMINATE, 2007")

# calculate pooled SD
dat$sigma2_placebo <- dat$sigma_placebo^2
dat$sigma2_treatment <- dat$sigma_treatment^2
dat$sigma2 <- dat$sigma2_placebo * (dat$n_placebo / dat$N) + dat$sigma2_treatment * (dat$n_treatment / dat$N)
dat$sigma <- sqrt(dat$sigma2)

# convert to SD units
results <- data.frame()
for (i in 1:nrow(dat)){
    rct <- get_rct_est(dat$mu_placebo[i], dat$n_placebo[i], dat$sigma_placebo[i], dat$mu_treatment[i], dat$n_treatment[i], dat$sigma_treatment[i])
    rct <- cbind(dat[i,], rct)
    results <- rbind(results, rct)
}

comb <- results %>% select(Target, Outcome) %>% unique
for (i in 1:nrow(comb)){
    # meta analysis of logSDR
    meta <- metagen(
        results %>% filter(Target == comb$Target[i] & Outcome == comb$Outcome[i]) %>% pull(est) %>% as.vector,
        results %>% filter(Target == comb$Target[i] & Outcome == comb$Outcome[i]) %>% pull(se) %>% as.vector,
        results %>% filter(Target == comb$Target[i] & Outcome == comb$Outcome[i]) %>% pull(Study) %>% as.vector,
        comb.fixed = TRUE, 
        comb.random = TRUE, 
        prediction=TRUE, 
        sm="SMD"
    )

    # forest plot
    png(paste0(comb$Target[i], "_", comb$Outcome[i], ".png"))
    forest(meta, layout = "JAMA")
    dev.off()
}