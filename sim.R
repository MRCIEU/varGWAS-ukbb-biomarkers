load("data/pheno.RData")
library('data.table')
library('dplyr')
library('optparse')
library('broom')
library("qqman")
source("funs.R")
set.seed(1234)

irnt <- function(pheno) {
	numPhenos = length(which(!is.na(pheno)))
	quantilePheno = (rank(pheno, na.last="keep", ties.method="random")-0.5)/numPhenos
	phenoIRNT = qnorm(quantilePheno)	
	return(phenoIRNT)
}

model <- function(data, out) {
  # prepare data
  data <- na.omit(data)
  data$snp <- get_simulated_genotypes(0.4, nrow(data))
  s <- "snp"
  data$xsq <- (data %>% pull(!!s))^2
  names(data) <- gsub("-", "_", names(data))
  out <- gsub("-", "_", out)

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
      Pvar=ftest$p.value[2]
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

opt <- data.frame(trait="body_mass_index.21001", stringsAsFactors=F)

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

# SD scale
dat[[opt$trait]] <- dat[[opt$trait]] / sd(dat[[opt$trait]])
dat$age_at_recruitment.21022.0.0 <- dat$age_at_recruitment.21022.0.0 / sd(dat$age_at_recruitment.21022.0.0)

# simulate exposure
n_obs <- 200
n_sim <- 5
af <- 0.4
results <- data.frame()
for (i in 1:n_sim){
    results <- rbind(results, model(dat, opt$trait))
}

# plot QQ
png(paste0("data/", opt$trait, "_phi_sim-null_qq.png"))
qq(results$P.r)
dev.off()