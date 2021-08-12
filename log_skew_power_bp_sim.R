library("broom")
library("metafor")
library("dplyr")
library("multcomp")
set.seed(123)

bp <- function(dat){
    dat$xsq <- dat$x^2
    fit1 <- lm(y ~ x, data=dat)
    dat$dsq <- resid(fit1)^2
    fit2 <- lm(dsq ~ x + xsq, data=dat)
    fitnull <- lm(dsq ~ 1, data=dat)
    return(tidy(anova(fitnull, fit2))$p.value[2])
}

n_obs <- 1000
n_sim <- 1000
m <- 1.44801
s <- 0.2

p_xmain <- rep(NA, n_sim)
p_umain <- rep(NA, n_sim)
p_int <- rep(NA, n_sim)
p_bp <- rep(NA, n_sim)
p_bpl <- rep(NA, n_sim)

for (i in 1:n_sim){
    x <- rbinom(n_obs, 1, 0.5)
    u <- rnorm(n_obs)
    y <- x*u*.12 + rlnorm(n=n_obs, meanlog=log(m^2 / sqrt(s^2 + m^2)), sdlog=sqrt(log(1 + (s^2 / m^2)))) # 80% pwr using B-P
    fit <- lm(y ~ x*u)
    p_xmain[i] <- tidy(fit)$p.value[2]
    p_umain[i] <- tidy(fit)$p.value[3]
    p_int[i] <- tidy(fit)$p.value[4]
    p_bp[i] <- bp(data.frame(x, u, y))
    p_bpl[i] <- bp(data.frame(x, u, y=log(y)))
}

binom.test(sum(p_xmain < 0.05), n_sim)
binom.test(sum(p_umain < 0.05), n_sim)
binom.test(sum(p_int < 0.05), n_sim)
binom.test(sum(p_bp < 0.05), n_sim)
binom.test(sum(p_bpl < 0.05), n_sim)