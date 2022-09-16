library("broom")
library("dplyr")

n_sim <- 200
n_obs <- 1000
set.seed(234)

results <- rep(NA, n_sim)
for (i in 1:n_sim){
    x <- rbinom(n_obs, 2, 0.25)
    u <- rbinom(n_obs, 2, 0.25)
    y <- x + u + x*u + rnorm(n_obs)
    r <- data.frame(x, u, y) %>% dplyr::group_by(x, u) %>% dplyr::summarize(mean=mean(y))
}