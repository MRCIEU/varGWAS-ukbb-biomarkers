load("data/pheno.RData")
library('optparse')
library('data.table')
library('dplyr')
library('broom')
library('ieugwasr')
library("lmtest")
library("sandwich")
source("funs.R")
set.seed(1234)
options(ieugwasr_api="http://web-dc1-bms-d0.infra.bris.ac.uk:5002/")

get_lad_resid <- function(X, y){
    X <- X %>% as.matrix
    y <- y %>% as.matrix
    fit <- cqrReg::qrfit(X=X, y=y, tau=.5, method="mm")
    b <- rbind(fit$b, fit$beta)
    X <- cbind(rep(1, nrow(X)), X)
    fitted <- X %*% b
    d <- y - fitted
    d <- abs(as.vector(d))
    return(d)
}

option_list = list(
  make_option(c("-t", "--trait"), type="character", default=NULL, help="Name of trait", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# export phenotypes for vGWAS
message(opt$trait)

# load linker
linker <- get_filtered_linker(drop_standard_excl=TRUE, drop_non_white_british=TRUE, drop_related=TRUE, application="15825")

# load covariates
covariates <- get_covariates()
covariates$chip <- as.numeric(as.factor(covariates$chip)) - 1
pc <- get_genetic_principal_components()

# merge data
dat <- merge(linker, covariates, "appieu")
dat <- merge(dat, pheno, by.x="app15825", by.y="eid")
dat <- merge(dat, pc, "appieu")

# SD scale env
for (e in env_exp){
  dat[[e]] <- dat[[e]] / sd(dat[[e]], na.rm=T)
}
dat[[opt$trait]] <- dat[[opt$trait]] / sd(dat[[opt$trait]], na.rm=T)

# read in GxE effects
gxe <- fread("data/gxe.txt")
gxe <- gxe %>% dplyr::filter(p.value < 5e-8)
gxe <- gxe %>% dplyr::filter(trait == !!opt$trait)
gxe <- cbind(gxe, str_split(gxe$term, ":", simplify=T) %>% as.data.frame %>% dplyr::rename(u="V1", x="V2"))
gxe <- cbind(gxe,  str_split(gxe$x, "_", simplify=T) %>% as.data.frame %>% dplyr::rename(chr="V1", pos="V2", oa="V3", ea="V4") %>% dplyr::mutate(pos=as.numeric(pos), chr=gsub("^chr", "", chr)))

# load dosages
for (i in 1:nrow(gxe)){
    dosage <- extract_variant_from_bgen(as.character(gxe$chr[i]), as.double(gxe$pos[i]), gxe$oa[i], gxe$ea[i])
    dat <- merge(dat, dosage, "appieu")
}

# test for interaction between each snp
results <- data.frame()
for (i in nrow(gxe)){
    X <- dat %>% 
        dplyr::select(age_at_recruitment.21022.0.0, sex.31.0.0, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, !!gxe$x[i], !!gxe$u[i], !!gxe$trait[i])
    X <- na.omit(X)
    y <- X %>% dplyr::select(!!gxe$trait[i])
    
    # test for effect on trait variance w/o adjusting for interaction effect
    X_temp <- X %>% dplyr::select(-!!gxe$trait[i], -!!gxe$u[i])
    d <- get_lad_resid(X_temp, y)
    X_temp$z <- as.factor(round(X_temp[[gxe$x[i]]]))
    X_temp <- X_temp %>% dplyr::select(-!!gxe$x[i])
    fit_unadj <- lm(d ~ ., data=X_temp)

    # test for effect on trait variance with adjusting for interaction effect
    X_temp <- X %>% dplyr::select(-!!gxe$trait[i])
    X_temp$xu <- X[[gxe$x[i]]] * X[[gxe$u[i]]]
    d <- get_lad_resid(X_temp, y)
    X_temp$z <- as.factor(round(X_temp[[gxe$x[i]]]))
    X_temp <- X_temp %>% dplyr::select(-!!gxe$x[i])
    fit_adj <- lm(d ~ ., data=X_temp)

    # compare model fits

    # save result
}