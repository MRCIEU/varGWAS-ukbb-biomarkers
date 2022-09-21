load("data/pheno.RData")

library("data.table")
library("ggplot2")
library("dplyr")
library("stringr")
library("ggpubr")
library("viridis")
library("lmtest")
library("sandwich")
library('forestplot')
library("viridis")
library("RColorBrewer")
library("grid")
library("broom")
library("gwasglue")
library("ieugwasr")
library("susieR")
library("robustbase")
source("funs.R")
options(ieugwasr_api="http://web-dc1-bms-d0.infra.bris.ac.uk:5002/")
set.seed(123)

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
    V1 <- as.data.frame(str_split(d$V1, "_", simplify=T), stringsAsFactors=F)
    names(V1) <- c("chr.1", "pos.1", "oa.1", "ea.1")
    V1$chr.1 <- gsub("chr", "", V1$chr.1)
    V2 <- as.data.frame(str_split(d$V2, "_", simplify=T), stringsAsFactors=F)
    names(V2) <- c("chr.2", "pos.2", "oa.2", "ea.2")
    V2$chr.2 <- gsub("chr", "", V2$chr.2)    
    d <- cbind(d, V1, V2)
    d$chr_pos.1 <- paste0(d$chr.1, ":", d$pos.1)
    d$chr_pos.2 <- paste0(d$chr.2, ":", d$pos.2)

    results <- data.frame()
    for (i in 1:nrow(d)){
        res <- all_interval_func(d$chr_pos.1[i], d$id[i])
        if (!is.null(res)){
            res$trait <- d$trait[i]
            results <- rbind(res, results)
        }
        res <- all_interval_func(d$chr_pos.2[i], d$id[i])
        if (!is.null(res)){
            res$trait <- d$trait[i]
            results <- rbind(res, results)
        }
    }

    return(results)
}

get_est <- function(trait, v1, v2, finemap){
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

    # SD normalise outcome
    dat[[trait]] <- dat[[trait]] / sd(dat[[trait]], na.rm=T)

    # subset finemapped variants for this trait
    finemap <- finemap %>% dplyr::filter(trait == !!trait)

    # load dosages
    snps1 <- stringr::str_split(v1, "_", simplify=T) %>% as.data.frame(., stringsAsFactors=F) %>% dplyr::rename(chr="V1", pos="V2", oa="V3", ea="V4") %>% dplyr::mutate(chr=gsub("chr", "", chr))
    snps2 <- stringr::str_split(v2, "_", simplify=T) %>% as.data.frame(., stringsAsFactors=F) %>% dplyr::rename(chr="V1", pos="V2", oa="V3", ea="V4") %>% dplyr::mutate(chr=gsub("chr", "", chr))
    snps3 <- finemap %>% dplyr::select(chr, position, oa, ea) %>% dplyr::rename(pos="position")
    snps <- unique(rbind(snps1, snps2, snps3))
    snps$key <- paste0("chr", snps$chr, "_", snps$pos, "_", snps$oa, "_", snps$ea)
    chr_pos.1 <- paste0(snps1$chr, ":", snps1$pos)
    chr_pos.2 <- paste0(snps2$chr, ":", snps2$pos)

    for (i in 1:nrow(snps)){
        dosage <- extract_variant_from_bgen(as.character(snps$chr[i]), as.double(snps$pos[i]), snps$oa[i], snps$ea[i])
        dat <- merge(dat, dosage, "appieu")
    }

    # test for interaction adjusting for finemapped variants
    adj1 <- finemap %>% dplyr::filter(chr_pos == !!chr_pos.1) %>% dplyr::select(chr, position, oa, ea) %>% dplyr::mutate(snp=paste0("chr", chr, "_", position, "_", oa, "_", ea)) %>% pull(snp)
    adj2 <- finemap %>% dplyr::filter(chr_pos == !!chr_pos.2) %>% dplyr::select(chr, position, oa, ea) %>% dplyr::mutate(snp=paste0("chr", chr, "_", position, "_", oa, "_", ea)) %>% pull(snp)
    adj <- unique(adj1, adj2)
    f <- paste0(trait, " ~ ", v1, " * ", v2, " + age_at_recruitment.21022.0.0 + sex.31.0.0 + ", paste0("PC", seq(1,10), collapse=" + "), " + ", paste0(adj, collapse=" + "))
    message(f)
    mod <- lm(f, data=dat)
    t <- coeftest(mod, vcov = vcovHC(mod, type = "HC0")) %>% tidy
    t <- t %>% dplyr::filter(grepl(":", t$term))
    t$trait <- trait
    t$formula <- f
    return(t)
}

# load gxe effects
additive <- get_dat("data/gxg.txt")

# append sensitivity P value
additive <- merge(additive, fread("data/gxg-log.txt") %>% dplyr::mutate(tt=paste0(trait, ":", term)) %>% dplyr::select(tt, p.value) %>% dplyr::rename(p_sens="p.value"), "tt")

# finemap gxg snps
additive.finemap_1 <- get_finemap(additive[1])
additive.finemap_1$target <- additive$V1[1]
additive.finemap_2 <- get_finemap(additive[2])
additive.finemap_2$target <- additive$V1[2]
additive.finemap_3 <- get_finemap(additive[3])
additive.finemap_3$target <- additive$V1[3]
additive.finemap_4 <- get_finemap(additive[4])
additive.finemap_4$target <- additive$V1[4]
additive.finemap_5 <- get_finemap(additive[5])
additive.finemap_5$target <- additive$V1[5]
additive.finemap_6 <- get_finemap(additive[6])
additive.finemap_6$target <- additive$V1[6]
additive.finemap_7 <- get_finemap(additive[7])
additive.finemap_7$target <- additive$V1[7]
additive.finemap_8 <- get_finemap(additive[8])
additive.finemap_8$target <- additive$V1[8]

additive.finemap <- rbind(
    additive.finemap_1,
    additive.finemap_2,
    additive.finemap_3,
    additive.finemap_4,
    additive.finemap_5,
    additive.finemap_6,
    additive.finemap_7,
    additive.finemap_8
)

# test for effect adjusting for finemapped variants
additive.results <- data.frame()
for (i in 1:nrow(additive)){
    est <- get_est(additive$trait[i], additive$V1[i], additive$V2[i], additive.finemap)
    additive.results <- rbind(additive.results, est)
}

write.table(additive.results, sep="\t", quote=F, row.names=F, file=paste0("data/gxg-add-finemap.txt"))
write.table(additive.finemap, sep="\t", quote=F, row.names=F, file=paste0("data/additive.finemap.txt"))