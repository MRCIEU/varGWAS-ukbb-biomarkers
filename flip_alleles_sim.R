library("broom")
library("dplyr")

n_obs <- 1000
n_sim <- 1000
set.seed(123)

bp <- function(x, y){
    x2 <- x^2
    fit1 <- lm(y ~ x)
    d <- resid(fit1)^2
    fit2 <- lm(d ~ x + x2)
    beta <- tidy(fit2)$estimate
    return(data.frame(x=beta[2], x2=beta[3]))
}

r1 <- data.frame()
r2 <- data.frame()
for (i in 1:n_sim){
    x <- rbinom(n_obs, 2, 1/3)
    u <- rnorm(n_obs)
    y <- x*u*.3 + rnorm(n_obs)
    r1 <- rbind(r1, bp(x,y))
    x <- dplyr::recode(x, `0` = 2L, `2` = 0L)
    r2 <- rbind(r2, bp(x,y))
}