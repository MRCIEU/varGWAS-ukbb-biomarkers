library("dplyr")
library("multcomp")
library("broom")
library("data.table")
library("stringr")
source("funs.R")
set.seed(1234)
options(ieugwasr_api="http://64.227.44.193:8006/")

bp <- function(dat, snp, outcome, covar){
    dat <- dat %>% dplyr::select(!!snp, !!outcome, !!covar) %>% tidyr::drop_na()
    dat[[outcome]] <- dat[[outcome]] / sd(dat[[outcome]])
    dat$x <- dat[[snp]]
    dat$xsq <- dat$x^2
    fit1 <- lm(paste0(outcome, " ~ x + ", paste0(covar, collapse= " + ")), data=dat)
    dat$dsq <- resid(fit1)^2
    fit2 <- lm(paste0("dsq ~ x + xsq + ", paste0(covar, collapse= " + ")), data=dat)
    fitnull <- lm(paste0("dsq ~ 1 + ", paste0(covar, collapse= " + ")), data=dat)
    varbeta1 <- glht(model=fit2, linfct=paste("x*1 + xsq*1 == 0"))
    varbeta2 <- glht(model=fit2, linfct=paste("x*2 + xsq*4 == 0"))
    res <- cbind(
        tidy(fit1) %>% dplyr::filter(term=="x") %>% dplyr::select("estimate", "std.error") %>% dplyr::rename(beta.estimate="estimate", beta.std.error="std.error"),
        tidy(varbeta1) %>% dplyr::select("estimate", "std.error") %>% dplyr::rename(varbeta1.estimate="estimate", varbeta1.std.error="std.error"),
        tidy(varbeta2) %>% dplyr::select("estimate", "std.error") %>% dplyr::rename(varbeta2.estimate="estimate", varbeta2.std.error="std.error"),
        n0=table(round(dat$x))[1],
        n1=table(round(dat$x))[2],
        n2=table(round(dat$x))[3],
        p=tidy(anova(fitnull, fit2))$p.value[2]
    )
    res$beta.lci <- res$beta.estimate - (res$beta.std.error * 1.96)
    res$beta.uci <- res$beta.estimate + (res$beta.std.error * 1.96)
    res$varbeta1.lci <- res$varbeta1.estimate - (res$varbeta1.std.error * 1.96)
    res$varbeta1.uci <- res$varbeta1.estimate + (res$varbeta1.std.error * 1.96)
    res$varbeta2.lci <- res$varbeta2.estimate - (res$varbeta2.std.error * 1.96)
    res$varbeta2.uci <- res$varbeta2.estimate + (res$varbeta2.std.error * 1.96)
    res$snp <- snp
    res$outcome <- outcome
    return(res)
}

# load phenotypes
load("data/pheno.RData")
linker <- get_filtered_linker(drop_standard_excl=TRUE, drop_non_white_british=TRUE, drop_related=TRUE, application="15825")
covariates <- get_covariates()
pc <- get_genetic_principal_components()
dat <- merge(linker, covariates, "appieu")
dat <- merge(dat, pheno, by.x="app15825", by.y="eid")
dat <- merge(dat, pc, "appieu")

# load clumped vQTLs
d <- fread("data/vqtls.txt")
d <- d %>% dplyr::filter(trait != "body_mass_index.21001.0.0")
d$key <- paste0("chr", d$chr, "_", d$pos, "_", d$oa, "_", d$ea)
d$chr_pos <- paste0(d$chr, ":", d$pos)

# load dosage
snps <- d %>% dplyr::select(chr, pos, oa, ea)
snps <- unique(snps)
for (i in 1:nrow(snps)){
    dosage <- tryCatch(
        expr = {
            extract_variant_from_bgen(as.character(snps$chr[i]), as.double(snps$pos[i]), snps$oa[i], snps$ea[i])
        },
        error = function(e){ 
            NULL
        }
    )

    if (is.null(dosage)){
        warning(paste0("skipping variant: ", snps$chr[i], ":", snps$pos[i]))
        next
    }
    dat <- merge(dat, dosage, "appieu")
}

results <- data.frame()
for (i in 1:nrow(d)){
    # main analysis
    res <- bp(dat, d$key[i], d$trait[i], c("age_at_recruitment.21022.0.0", "sex.31.0.0", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"))

    # log-scale analysis
    dat[[paste0(d$trait[i], "_log")]] <- log(dat[[d$trait[i]]])
    res_log <- bp(dat, d$key[i], paste0(d$trait[i], "_log"), c("age_at_recruitment.21022.0.0", "sex.31.0.0", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"))
    names(res_log) <- paste0(names(res_log), ".log")

    results <- rbind(results, cbind(res, res_log))
}

# select fields for paper
results <- results %>% dplyr::select("snp", "outcome", "beta.estimate", "beta.lci", "beta.uci", "varbeta1.estimate", "varbeta1.lci", "varbeta1.uci", "varbeta2.estimate", "varbeta2.lci", "varbeta2.uci", "p", "beta.estimate.log", "beta.lci.log", "beta.uci.log", "varbeta1.estimate.log", "varbeta1.lci.log", "varbeta1.uci.log", "varbeta2.estimate.log", "varbeta2.lci.log", "varbeta2.uci.log", "p.log")
results$key <- stringr::str_split(results$snp, "_", simplify=T) %>% as.data.frame %>% dplyr::mutate(V1=gsub("chr", "", V1)) %>% dplyr::mutate(key=paste0(V1, ":",V2)) %>% dplyr::pull(key)

# get data on nearest gene
ng <- fread("data/nearest.txt")
ng <- unique(ng)
ng$key <- paste0(ng$V1,":",ng$V2)
ng$V1 <- NULL
ng$V2 <- NULL
names(ng)[1] <- "gene"
ng <- ng %>%
  group_by_at(vars(key)) %>%
  summarize(gene = toString(gene)) %>%
  ungroup()
ng$gene <- gsub(", ", "|", ng$gene)

# results with nearest gene
all <- merge(results, ng, "key", all.x=T)
all$key <- NULL

# tidy outcome name
all$outcome <- sapply(all$outcome, function(x) biomarkers_abr[x==biomarkers], simplify=T)

# write to table
write.csv(results, file="Table S1.csv", quote=F, row.names=F)