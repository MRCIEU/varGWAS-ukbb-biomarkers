load("data/pheno.RData")
library('data.table')
library('dplyr')
library('optparse')
library('broom')
library("qqman")
library("genpwr")
source("funs.R")
set.seed(1234)

option_list = list(
  make_option(c("-t", "--trait"), type="character", default=NULL, help="Name of trait", metavar="character"),
  make_option(c("-n", "--norm"), type="character", default=NULL, help="Normlisation strategy", metavar="character"),
  make_option(c("-m", "--mean"), type="character", default=NULL, help="Mean effect", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

irnt <- function(pheno) {
	numPhenos = length(which(!is.na(pheno)))
	quantilePheno = (rank(pheno, na.last="keep", ties.method="random")-0.5)/numPhenos
	phenoIRNT = qnorm(quantilePheno)	
	return(phenoIRNT)
}

model <- function(data, out, n_obs, af, delta) {
  # prepare data
  data <- data[sample(nrow(data), n_obs, replace=F), ]
  data$snp <- get_simulated_genotypes(af, nrow(data))
  s <- "snp"
  data$xsq <- (data %>% pull(!!s))^2
  names(data) <- gsub("-", "_", names(data))
  out <- gsub("-", "_", out)

  # set mean effect of snp
  data[[out]] <- data[[s]] * delta + data[[out]]

  # first-stage model
  f <- paste0(out, " ~ ", s, " + sex.31.0.0 + age_at_recruitment.21022.0.0 +", paste0("PC", seq(1, 10), collapse="+"))
  fit1 <- lm(as.formula(f), data=data)

  # second-stage model
  data$d <- resid(fit1)^2
  f <- paste0("d ~ ", s, " + xsq")
  fit2 <- lm(as.formula(f), data=data)

  # F-test
  fit0 <- lm(d ~ 1, data=data)
  ftest <- tidy(anova(fit0, fit2))
  fit2t <- tidy(fit2)

  return(data.frame(
      SNP=s,
      BETA_x=fit2t$estimate[2],
      BETA_xq=fit2t$estimate[3], 
      SE_x=fit2t$std.error[2],
      SE_xq=fit2t$std.error[3],
      Pvar=ftest$p.value[2],
      Pmu=tidy(fit1)$p.value[2]
    )
  )
}

#' Function to simulate genotypes in HWE
#' @param q Recessive/alternative allele frequency
#' @param n_obs Number of observations to return
get_simulated_genotypes <- function(q, n_obs) {
  p <- 1 - q
  x <- sample(c(0, 1, 2), n_obs, prob = c(p^2, 2 * p * q, q^2), replace = T)
  return(x)
}

n_sim <- 1000
n_obs <- 100000
af <- 0.4

# load linker
linker <- get_filtered_linker(drop_standard_excl=TRUE, drop_non_white_british=TRUE, drop_related=TRUE, application="15825")

# load covariates
covariates <- get_covariates()
pc <- get_genetic_principal_components()

# merge data
dat <- merge(linker, covariates, "appieu")
dat <- merge(dat, pheno, by.x="app15825", by.y="eid")
dat <- merge(dat, pc, "appieu")

# select fields for GWAS
dat <- dat[,c("appieu", "sex.31.0.0", "age_at_recruitment.21022.0.0", opt$trait, paste0("PC", seq(1, 10))), with=F]

# drop missing values
dat <- dat[complete.cases(dat), ]

# standardise
dat$age_at_recruitment.21022.0.0 <- dat$age_at_recruitment.21022.0.0 / sd(dat$age_at_recruitment.21022.0.0)

if (opt$norm == "sd"){
    dat[[opt$trait]] <- dat[[opt$trait]] / sd(dat[[opt$trait]])
} else if (opt$norm == "irnt"){
    dat[[opt$trait]] <- irnt(dat[[opt$trait]])
} else {
    stop(paste0("Norm strategy not known:", opt$norm))
}

if (opt$mean == "True"){
    # main effect size of X on Y detectable with 80% power
    delta <- as.numeric(genpwr.calc(calc = "es", model = "linear", ge.interaction = NULL, N = n_obs, k = NULL, MAF = af, Power = 0.8, Alpha = 0.05, sd_y = sd(dat[[opt$trait]]), True.Model = "Additive", Test.Model = "Additive")$ES_at_Alpha_0.05)
    message(paste0("setting main effect to ", delta))
}

# simulate exposure with no effect on outcome & test for variance effect
results <- data.frame()
for (i in 1:n_sim){
    message(paste0("n_sim: ", i))
    results <- rbind(results, model(dat, opt$trait, n_obs, af, delta))
}

# plot QQ
png(paste0("data/", opt$trait, "_sim_", opt$norm, "_", opt$mean,"_phi_qq.png"))
qq(results$Pvar)
dev.off()

png(paste0("data/", opt$trait, "_sim_", opt$norm, "_", opt$mean,"_mu_qq.png"))
qq(results$Pmu)
dev.off()

# T1E
write.csv(tidy(binom.test(sum(results$Pvar < 0.05), n_sim)), row.names=F, file=paste0("data/", opt$trait, "_sim_", opt$norm, "_", opt$mean,"_phi_t1e.csv"))
write.csv(tidy(binom.test(sum(results$Pmu < 0.05), n_sim)), row.names=F, file=paste0("data/", opt$trait, "_sim_", opt$norm, "_", opt$mean,"_mu_t1e.csv"))