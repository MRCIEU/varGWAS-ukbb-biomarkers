library("metafor")
library("readxl")
library('meta')
library("dplyr")
library("data.table")
library("ggplot2")
set.seed(1234)

### meta analysis with log CVR
# NB note that escalc ignores the rhos from the sampling variance, hence assumes normality of data
# outputs data.frame with columns of logCVR, SE, lower and upper CI, test statistic and pvalue
# Taken from https://github.com/harrietlmills/DetectingDifferencesInVariance/blob/e8cbb3426f5475c027b3e6e08fc19ee811b3b759/MetaAnalysis.R#L236
MA_analysis_logCVR <- function(MA_dataset){
  # MA_dataset is a dataframe with columns: "Study"     "Int_Mean"  "Int_SD"    "Int_N"     "Con_Mean"  "Con_SD"    "Con_N" 
  
  # calculate CVR
  rdat <- escalc(measure = "CVR", 
                 m1i = MA_dataset$Int_Mean, n1i = MA_dataset$Int_N, sd1i = MA_dataset$Int_SD, 
                 m2i = MA_dataset$Con_Mean, n2i = MA_dataset$Con_N, sd2i = MA_dataset$Con_SD)
  
  # calculate confidence intervals, etc
  srdat <- summary(rdat, digits = 2)
  # yi-observed outcome, vi-estimated sampling vars, sei-standard errors of observed outcomes, zi=test statistics, ci.lb and ci.ub=conf intervals,
  # if wanted, can specify trans = exp, which means the observed outcome and upper and lower bounds are transformed
  
  # fit random-effects models
  m1    <- rma(yi = yi, vi = vi, data = rdat, method = "REML",
               slab = MA_dataset$Study, weighted = TRUE)
  sum_m1<-summary(m1)
  #coef(summary(m1))
  
  
  MA_table <- cbind(srdat$yi, srdat$sei, srdat$ci.lb, srdat$ci.ub, srdat$zi, NA) #srdat$vi)
  rownames(MA_table) <- MA_dataset$Study
  colnames(MA_table) <- paste0("logCVR_", c("est", "SE", "Lower", "Upper", "zval", "pval"))
  
  random <- c(sum_m1$beta, sum_m1$se, sum_m1$ci.lb, sum_m1$ci.ub, sum_m1$zval, sum_m1$pval)
  MA_table <- rbind(MA_table, logCVR_random=random)
  
  res <- data.frame(MA_table)
  res$Study <- rownames(res)
  return(res)
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

# restrict to required fields
rct <- dat %>% select(Outcome, Study, mu_treatment, sigma_treatment, n_treatment, mu_placebo, sigma_placebo, n_placebo) %>% 
    rename(Int_Mean=mu_treatment, Int_SD=sigma_treatment, Int_N=n_treatment, 
        Con_Mean=mu_placebo, Con_SD=sigma_placebo, Con_N=n_placebo)

# MR
dat <- fread("iv_snp_lipid_variance.txt")
dat$outcome <- gsub("ldl_direct.30780.0.0", "LDL", dat$outcome)
dat$outcome <- gsub("hdl_cholesterol.30760.0.0", "HDL", dat$outcome)
mr <- dat %>% select(outcome, snp, estimate.2, sd.estimate.2, n2, estimate.0, sd.estimate.0, n0) %>% 
    rename(Outcome=outcome, Study=snp, Int_Mean=estimate.2, Int_SD=sd.estimate.2, Int_N=n2, 
        Con_Mean=estimate.0, Con_SD=sd.estimate.0, Con_N=n0)

# meta-analsyis of RCT
ldl.rct <- MA_analysis_logCVR(rct %>% filter(Outcome == "LDL"))
hdl.rct <- MA_analysis_logCVR(rct %>% filter(Outcome == "HDL"))

# meta-analsyis of MR
ldl.mr <- MA_analysis_logCVR(mr %>% filter(Outcome == "LDL"))
hdl.mr <- MA_analysis_logCVR(mr %>% filter(Outcome == "HDL"))

# combined meta-analysis
ldl.both <- MA_analysis_logCVR(rbind(rct %>% filter(Outcome == "LDL"), mr %>% filter(Outcome == "LDL")))
hdl.both <- MA_analysis_logCVR(rbind(rct %>% filter(Outcome == "HDL"), mr %>% filter(Outcome == "HDL")))

# forest plot
get_plot <- function(d){
    d$logCVR_est <- exp(d$logCVR_est)
    d$logCVR_Lower <- exp(d$logCVR_Lower)
    d$logCVR_Upper <- exp(d$logCVR_Upper)
    p <- ggplot(d, aes(x=Study, y=logCVR_est, ymin=logCVR_Lower, ymax=logCVR_Upper)) +
        coord_flip() +
        scale_colour_brewer(palette = "Set1") +
        geom_point() +
        geom_errorbar(width=.05) +
        geom_hline(yintercept = 1, linetype = "dashed", color = "grey") +
        theme_classic()
    return(p)
}