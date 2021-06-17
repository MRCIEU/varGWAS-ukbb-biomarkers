library("jlst")
library("data.table")
library("pwr")
library("broom")
library("qqman")
source("funs.R")
set.seed(123)

#' Estimate treatment effect size a categorical exposure
#' @param n_obs Number of samples in analysis
#' @param size Number of categorical levels - 1
#' @param treatment1_p The probability of recieving treatment
#' @param sd The SD out of the normal outcome
#' @param alpha P value threshold
#' @param power Power to detect effect
#' @return Effect size
get_cat_delta <- function(n_obs, size, treatment1_p, sd, alpha, power){
    p <- power.t.test(n=n_obs * (size * treatment1_p), delta=NULL, sd=sd, sig.level=alpha, power=power, type = c("two.sample"), alternative = c("two.sided"))
    return(p$delta)
}

n_sim <- 1000
n_obs <- 10000

# read in outcome
dat <- fread("data/aspartate_aminotransferase.30650.0.0.txt")

# calculate main effect size for categorical exposure to detect at specified power
b <- get_cat_delta(n_obs, 2, 0.5, 1, 0.05, 0.8)
b <- 0

# simulate SNP effect
p_var <- rep(NA, n_sim)
p_int <- rep(NA, n_sim)
for (i in 1:n_sim){
    x <- rbinom(n=n_obs, size=2, prob=.5)
    u <- rbinom(n=n_obs, size=2, prob=.5)
    y <- x*b + sample(dat$aspartate_aminotransferase.30650.0.0, n_obs)
    t <- vartest(y, x, type = 1, x.sq = T)
    p_var[i] <- t$test$P
    p_int[i] <- (tidy(lm(y~x*u)) %>% pull(p.value))[4]
}

# estimate power
tidy(binom.test(sum(p_var < 0.05), n_sim))
tidy(binom.test(sum(p_int < 0.05), n_sim))

# plot qq
qq(p_var)
qq(p_int)