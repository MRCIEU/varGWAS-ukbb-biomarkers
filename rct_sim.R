library("broom")
library("metafor")
library("dplyr")
set.seed(123)

# LogCVR = t_cv / c_cv
# CVR < 1 = SD smaller in treatment arm
# CVR > 1 = SD larger in treatment arm

get_logcvr <- function(y0, y1){
    # Taken from https://github.com/harrietlmills/DetectingDifferencesInVariance/blob/e8cbb3426f5475c027b3e6e08fc19ee811b3b759/AnalyseIndividualTrials.R#L253
    ### log of ratio of CoVs (logCVR)
    rdat_logCVR <- escalc(measure = "CVR", #m1=exp, m2=con
                        m1i = mean(y1), n1i = length(y1), sd1i = sd(y1), 
                        m2i = mean(y0), n2i = length(y0), sd2i = sd(y0))

    # calculate confidence intervals, etc
    srdat_logCVR <- summary(rdat_logCVR, digits = 4) 

    # calculate test statistic and pvalue
    logCVR_pvalue <- 2*(1-pnorm(abs(srdat_logCVR$zi))) 

    return(data.frame(
        estimate=srdat_logCVR$yi,
        std.error=srdat_logCVR$sei,
        t.statistic=srdat_logCVR$zi,
        p.value=logCVR_pvalue
    ))
}

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
  
  return(data.frame(MA_table))
}

n_obs <- 1000
n_sim <- 1000
m <- 1.44801 # mean of HDL-c in UKBB
s <- 0.382304 # sd of HDL-c in UKBB

# RCT effect on cholesterol conc
p_xmain <- rep(NA, n_sim)
p_umain <- rep(NA, n_sim)
p_int <- rep(NA, n_sim)
results <- data.frame()
ma <- data.frame()
for (i in 1:n_sim){
    x <- rbinom(n_obs, 1, 0.5)
    u <- rnorm(n_obs)
    y <- x*.175 + x*u*.175+ rlnorm(n=n_obs, meanlog=log(m^2 / sqrt(s^2 + m^2)), sdlog=sqrt(log(1 + (s^2 / m^2)))) # effects set to have 80% pwr
    fit <- lm(y ~ x*u)
    p_xmain[i] <- tidy(fit)$p.value[2]
    p_umain[i] <- tidy(fit)$p.value[3]
    p_int[i] <- tidy(fit)$p.value[4]
    results <- rbind(results, get_logcvr(y[x==0], y[x==1]))
    ma <- rbind(ma, data.frame(Study=i, Int_Mean=mean(y[x==1]), Int_SD=sd(y[x==1]), Int_N=length(y[x==1]), Con_Mean=mean(y[x==0]), Con_SD=sd(y[x==0]), Con_N=length(y[x==0])))
}

binom.test(sum(p_xmain < 0.05), n_sim)
binom.test(sum(p_umain < 0.05), n_sim)
binom.test(sum(p_int < 0.05), n_sim)
MA_analysis_logCVR(ma)

# GxE effect
q <- 1 / 3
p_xmain <- rep(NA, n_sim)
p_umain <- rep(NA, n_sim)
p_int <- rep(NA, n_sim)
for (i in 1:n_sim){
    x <- rbinom(n_obs, 2, q)
    u <- rnorm(n_obs)
    y <- x*.1298 + x*u*.1345 + rnorm(n_obs)
    fit <- lm(y ~ x*u)
    p_xmain[i] <- tidy(fit)$p.value[2]
    p_umain[i] <- tidy(fit)$p.value[3]
    p_int[i] <- tidy(fit)$p.value[4]
} 

# power
binom.test(sum(p_xmain < 0.05), n_sim)
binom.test(sum(p_umain < 0.05), n_sim)
binom.test(sum(p_int < 0.05), n_sim)
