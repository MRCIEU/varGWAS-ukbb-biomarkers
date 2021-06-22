library("jlst")
library("data.table")
library("pwr")
library("broom")
library('optparse')
library("qqman")
source("funs.R")
set.seed(123)

option_list = list(
  make_option(c("-t", "--trait"), type="character", default=NULL, help="Name of trait", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#' Function to simulate genotypes in HWE
#' @param q Recessive/alternative allele frequency
#' @param n_obs Number of observations to return
get_simulated_genotypes <- function(q, n_obs){
  p <- 1 - q
  x <- sample(c(0, 1, 2), n_obs, prob=c(p^2, 2 * p * q, q^2), replace=T)
  return(x)
}

# read in outcome
dat <- fread(paste0("data/", opt$trait, ".txt")) %>% 
    pull(opt$trait)

n_sim <- 1000
n_obs <- length(dat)

# simulate SNP effect
results <- data.frame()
for (q in c(0.05, 0.1, 0.15, 0.2, 0.25)){
    for (i in 1:n_sim){
        x <- get_simulated_genotypes(q, n_obs)
        y <- x*0 + dat
        t <- vartest(y, x, type = 1, x.sq = T)
        results <- rbind(results, data.frame(q, p=t$test$P))
    }
}

# estimates T1E for each value of q
results <- results %>%
    group_by(q) %>%
    summarize(tidy(binom.test(sum(p < 0.05), n_sim)))

# write out results
write.table(results, sep="\t", quote=F, row.names=F, file=paste0("data/", opt$trait, ".null.txt"))