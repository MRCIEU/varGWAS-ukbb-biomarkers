set.seed(123)

# test effect of RCT mean increase on SD

# params from Lincoff et al, 2017
n_p <- 6054
n_t <- 6038
n <- 6054 + 6038
x <- sample(c(0, 1), n, c(6054/n, 6038/n), replace=T)
y <- 45.6 + x*58.5 + rnorm(n, sd=12.3)
sd(y[x==0])
sd(y[x==1])