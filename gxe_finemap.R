load("data/pheno.RData")

library("data.table")
library("ggplot2")
library("dplyr")
library("stringr")
library("ggpubr")
library("viridis")
library("lmtest")
library("sandwich")
library("viridis")
library("RColorBrewer")
library("grid")
library("broom")
library("optparse")
library("gwasglue")
library("ieugwasr")
library("susieR")
library("robustbase")
source("funs.R")
options(ieugwasr_api="http://64.227.44.193:8006/")
set.seed(123)

option_list = list(
  make_option(c("-t", "--trait"), type="character", default=NULL, help="Name of trait", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

get_dat <- function(file){
    # read in gxe results
    d <- fread(file)

    # drop BMI
    d <- d %>% dplyr::filter(trait != "body_mass_index.21001.0.0")

    # filter SNPs to show
    d <- d %>% dplyr::filter(p.value < 5e-8)

    # merge
    d <- cbind(d, as.data.frame(stringr::str_split(d$term, ":", simplify=T), stringsAsFactors=F))
    
    # add key
    d$tt <- paste0(d$trait, ":", d$term)

    return(d)
}

get_finemap <- function(d){
    # add opengwas ID
    d$id <- paste0("ukb-d-", str_split(d$trait, "\\.", simplify=T)[,2], "_irnt")

    # extract chr-pos
    V2 <- as.data.frame(str_split(d$V2, "_", simplify=T), stringsAsFactors=F)
    names(V2) <- c("chr.1", "pos.1", "oa.1", "ea.1")
    V2$chr.1 <- gsub("chr", "", V2$chr.1)
    d <- cbind(d, V2)
    d$chr_pos.1 <- paste0(d$chr.1, ":", d$pos.1)

    results <- data.frame()
    for (i in 1:nrow(d)){
        res <- finemap_func(d$chr_pos.1[i], d$id[i])
        if (!is.null(res)){
            res$trait <- d$trait[i]
            res$target <- d$V2[i]
            res$modifier <- d$V1[i]
            results <- rbind(res, results)
        }
    }

    return(results)
}

get_est <- function(trait, v1, v2, finemap, pheno){
    # subset finemapped variants for this trait
    finemap <- finemap %>% dplyr::filter(trait == !!trait)
    finemap$key <- paste0("chr", finemap$chr, "_", finemap$pos, "_", finemap$oa, "_", finemap$ea)

    # load dosages
    snps2 <- stringr::str_split(v2, "_", simplify=T) %>% as.data.frame(., stringsAsFactors=F) %>% dplyr::rename(chr="V1", pos="V2", oa="V3", ea="V4") %>% dplyr::mutate(chr=gsub("chr", "", chr))
    snps3 <- finemap %>% dplyr::select(chr, position, oa, ea) %>% dplyr::rename(pos="position")
    snps <- unique(rbind(snps2, snps3))
    snps$key <- paste0("chr", snps$chr, "_", snps$pos, "_", snps$oa, "_", snps$ea)
    chr_pos.2 <- paste0(snps2$chr, ":", snps2$pos)

    skiped <- data.frame()
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
            skiped <- rbind(skiped, data.frame(key=paste0("chr", snps$chr[i], "_", snps$pos[i], "_", snps$oa[i], "_", snps$ea[i])))
            next
        }
        pheno <- merge(pheno, dosage, "appieu")
    }

    # test for interaction adjusting for finemapped variants
    adj <- finemap %>% dplyr::filter(chr_pos == !!chr_pos.2) %>% dplyr::select(chr, position, oa, ea) %>% dplyr::mutate(snp=paste0("chr", chr, "_", position, "_", oa, "_", ea)) %>% pull(snp) %>% unique
    adj <- adj[!adj %in% skiped$key]
    if (length(adj) > 0){
        f <- paste0(trait, " ~ ", v1, " * ", v2, " + age_at_recruitment.21022.0.0 + sex.31.0.0 + ", paste0("PC", seq(1,10), collapse=" + "), " + ", paste0(adj, collapse=" + "))
    } else {
        f <- paste0(trait, " ~ ", v1, " * ", v2, " + age_at_recruitment.21022.0.0 + sex.31.0.0 + ", paste0("PC", seq(1,10), collapse=" + "))
    }
    
    message(f)
    mod <- lm(f, data=pheno)
    t <- coeftest(mod, vcov = vcovHC(mod, type = "HC0")) %>% tidy
    t <- t %>% dplyr::filter(grepl(":", t$term))
    t$trait <- trait
    t$formula <- f
    return(t)
}

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

# load gxe effects
additive <- get_dat("data/gxe.txt")

# append sensitivity P value
additive <- merge(additive, fread("data/gxe-log.txt") %>% dplyr::mutate(tt=paste0(trait, ":", term)) %>% dplyr::select(tt, p.value) %>% dplyr::rename(p_sens="p.value"), "tt")

# subet for this trait
additive <- additive %>% dplyr::filter(trait == !!opt$trait)

# stop if no GxE effects
stopifnot(nrow(additive) > 0)

# finemap gxg snps
additive.finemap <- get_finemap(additive)

# test for effect adjusting for finemapped variants
additive.results <- data.frame()
for (i in 1:nrow(additive)){
    est <- get_est(additive$trait[i], additive$V1[i], additive$V2[i], additive.finemap, dat)
    additive.results <- rbind(additive.results, est)
}

write.table(additive.results, sep="\t", quote=F, row.names=F, file=paste0("data/", opt$t,"_gxe-add-finemap.txt"))
write.table(additive.finemap, sep="\t", quote=F, row.names=F, file=paste0("data/", opt$t,"_gxe-additive.finemap.txt"))