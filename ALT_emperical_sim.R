library("dplyr")
library("broom")
library("pwr")
library("jlst")
library("tidyr")
library("ggpubr")
library("data.table")
set.seed(123)

n_sim <- 1000
n_obs <- 1000

# read in emperical dist
d <- fread("data/alanine_aminotransferase.30620.0.0.txt")

#' Function to estimate the power of an MC experiment
#' @param results Dataframe containing: P value and analysis group(s)
#' @param field Name of P value to summarise
#' @param n_sim Number of simulations performed
#' @param grp_name A vector of fields for use in grouping the analysis
#' @param alpha P threshold
calc_power <- function(results, field, n_sim, grp_name, alpha=0.05){
    # threshold P value
    results <- results %>%
        mutate(pos = as.numeric(get(field) < alpha))
    
    # count positives and estimate power with 95% CI
    h1 <- results %>%
        select(pos, all_of(grp_name)) %>%
        drop_na() %>%
        group_by_at(vars(all_of(grp_name))) %>%
        summarise(h1 = sum(pos)) %>%
        rowwise() %>%
        mutate(est_power = tidy(binom.test(h1, n_sim))$estimate) %>% 
        mutate(est_power_low = tidy(binom.test(h1, n_sim))$conf.low) %>%
        mutate(est_power_high = tidy(binom.test(h1, n_sim))$conf.high)

    # factorise grouping variables
    h1[grp_name] <- lapply(h1[grp_name], factor)

    return(h1)
}

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

#' Estimate treatment effect size given a continuous exposure
#' @param n_obs Number of samples in analysis
#' @param alpha P value threshold
#' @param power Power to detect effect
#' @return Effect size
get_cont_delta <- function(n_obs, alpha, power){
    # u = number of terms in the model (excluding intercept)
    u = 1

    # v = error degrees of freedom
    v = n_obs - u - 1

    # f2 = Cohen f2 (variance explained by the model)
    p <- pwr.f2.test(u = u, v = v, sig.level = alpha, power = power)

    return(sqrt(p$f2))
}

ols <- function(x, y){
    fit <- lm(y ~ x)
    fit <- tidy(fit)
    return(fit$p.value[2])
}

# QQ plots
qqgplot <- function(data, test, exposure, ci = 0.95) {
    temp <- data.frame()

    # straify by exposure variable
    data <- data %>%
        filter(exposure == !!exposure)

    # estimate power/fpr
    for (effect in c("TP", "TN")){
        for (dist in c("Emperical")){
            p <- data %>%
                filter(dist == !!dist & effect == !!effect) %>%
                pull(!!test)
            n  <- length(p)
            temp <- rbind(temp, data.frame(
                dist,
                effect,
                observed = -log10(sort(p)),
                expected = -log10(ppoints(n)),
                clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
                cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
            ))
        }
    }
    temp$dist <- factor(temp$dist, levels = c("Normal", "T", "Lognormal", "Mixed Normal"))
    temp$effect <- factor(temp$effect, levels = c("TP", "TN"))

    log10Pe <- expression(paste("Expected -log"[10], plain(P)))
    log10Po <- expression(paste("Observed -log"[10], plain(P)))

    pl <- ggplot(temp) +
        geom_point(aes(expected, observed), shape = 1, size = 1) +
        geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
        geom_line(aes(expected, cupper), linetype = 2) +
        geom_line(aes(expected, clower), linetype = 2) +
        xlab(log10Pe) +
        ylab(log10Po) +
        facet_grid(vars(effect), vars(dist), scales = "free_y") +
        theme_minimal()

    return(pl)
}

results <- data.frame()
for (effect in c("TP", "TN")){
    for (exposure in c("Binary", "Categorical", "Continuous")){
        for (i in 1:n_sim){
            if (exposure == "Binary"){
                if (effect == "TP"){
                    # calculate main effect size for categorical exposure to detect at specified power
                    b <- get_cat_delta(n_obs, 1, 0.5, 1, 0.05, 0.8)
                } else {
                    b <- 0
                }
                x <- rbinom(n_obs, 1, 0.5)
            } else if (exposure == "Categorical") {
                if (effect == "TP"){
                    # calculate main effect size for categorical exposure to detect at specified power
                    b <- get_cat_delta(n_obs, 2, 0.5, 1, 0.05, 0.8)
                } else {
                    b <- 0
                }
                x <- rbinom(n_obs, 2, 0.5)
            } else if (exposure == "Continuous") {
                if (effect == "TP"){
                    # calculate main effect size for continuous exposure to detect at specified power
                    b <- get_cont_delta(n_obs, 0.05, 0.8)
                } else {
                    b <- 0
                }
                x <- rnorm(n_obs)
            }
            
            y <- x * b + sample(d$alanine_aminotransferase.30620.0.0, n_obs, replace=T)

            # test for effect using battery of tests
            res <- data.frame(
                effect,
                exposure,
                dist="Emperical",
                b,
                ols_p=ols(x, y),
                bp_p=vartest(y, x, type=1)$test$P,
                bf_p=vartest(y, x, type=2)$test$P,
                jlsp_p=jlsp(y, x, var.type=2)$location_scale_test$P,
                jlssc_p=jlssc(y, x, type=2)$P
            )
            results <- rbind(results, res)
        }
    }
}

# estimate false-positive rate and power for each analysis
ols_est <- calc_power(results, "ols_p", n_sim, c("effect", "exposure", "dist"), 0.05)
write.table(ols_est, quote=F, row.names=F, sep="\t", file="ols_power.txt")

bp_est <- calc_power(results, "bp_p", n_sim, c("effect", "exposure", "dist"), 0.05)
write.table(bf_est, quote=F, row.names=F, sep="\t", file="bp_power.txt")

bf_est <- calc_power(results, "bf_p", n_sim, c("effect", "exposure", "dist"), 0.05)
write.table(bf_est, quote=F, row.names=F, sep="\t", file="bf_power.txt")

# binary exposure
p1 <- qqgplot(results, "ols_p", "Binary")
p2 <- qqgplot(results, "bf_p", "Binary")
p3 <- qqgplot(results, "bp_p", "Binary")

bin_p <- ggarrange(p1, p2, p3, labels = c("A", "B", "C"), ncol = 2, nrow = 2)
pdf("bin_p.pdf")
print(bin_p)
dev.off()

# categorical exposure
p1 <- qqgplot(results, "ols_p", "Categorical")
p2 <- qqgplot(results, "bf_p", "Categorical")
p2 <- qqgplot(results, "bp_p", "Categorical")

cat_p <- ggarrange(p1, p2, p3, labels = c("A", "B", "C"), ncol = 2, nrow = 2)
pdf("cat_p.pdf")
print(cat_p)
dev.off()

# continuous exposure
p1 <- qqgplot(results, "ols_p", "Continuous")
p2 <- qqgplot(results, "bf_p", "Continuous")
p2 <- qqgplot(results, "bp_p", "Continuous")

cont_p <- ggarrange(p1, p2, p3, labels = c("A", "B", "C"), ncol = 2, nrow = 2)
pdf("cont_p.pdf")
print(cont_p)
dev.off()