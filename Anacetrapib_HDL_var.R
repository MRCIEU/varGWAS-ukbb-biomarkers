# Taken from https://stats.stackexchange.com/questions/256456/how-to-calculate-mean-and-standard-deviation-from-median-and-quartiles
# Taken from https://www.ahajournals.org/doi/10.1161/JAHA.120.018136
get_sd <- function(q1, q3, n){
    s <- (q3 - q1) / (2 * (qnorm((0.75 * n - 0.125) / (n + 0.25))))
    return(s)
}

# placebo
placebo_n <- 289
placebo_mean <- (38+45+51) / 3
placebo_sd <- get_sd(38, 51, 289)

# treatment
treatment_n <- 285
treatment_mean <- (82+100+117) / 3
treatment_sd <- get_sd(82, 117, 285)

# test for equal variance
# df=1, T=244.620276094346, P=0

### log of ratio of variability (logVR)
library(metafor)

# Taken from https://github.com/harrietlmills/DetectingDifferencesInVariance/blob/master/AnalyseIndividualTrials.R
rdat_logVR <- escalc(measure = "VR", #m1=exp, m2=con
                     m1i = placebo_mean, n1i = placebo_n, sd1i = placebo_sd, 
                     m2i = treatment_mean, n2i = treatment_n, sd2i = treatment_sd)

# calculate confidence intervals, etc
srdat_logVR <- summary(rdat_logVR, digits = 2)

# calculate test statistic and pvalue
logVR_pvalue <-  2*(1-pnorm(abs(srdat_logVR$zi)))

est <- rdat_logVR$yi 
se <-  srdat_logVR$sei 
test <- srdat_logVR$zi
p <- logVR_pvalue