library("ieugwasr")
library("dplyr")
library("multcomp")
library("broom")
library("tidyverse")
library("data.table")
library("IRanges")
library("stringr")
library("GenomicRanges")
library("TwoSampleMR")
library("susieR")
library("gwasglue")
source("funs.R")
set.seed(1234)
options(ieugwasr_api="http://64.227.44.193:8006/")

bp <- function(dat, snp, outcome, covar){
    dat <- dat %>% dplyr::select(!!snp, !!outcome, !!covar) %>% tidyr::drop_na()
    dat$x <- dat[[snp]]
    dat$xsq <- dat$x^2
    fit1 <- lm(paste0(outcome, " ~ x + ", paste0(covar, collapse= " + ")), data=dat)
    dat$dsq <- resid(fit1)^2
    fit2 <- lm(paste0("dsq ~ x + xsq + ", paste0(covar, collapse= " + ")), data=dat)
    fitnull <- lm(paste0("dsq ~ 1 + ", paste0(covar, collapse= " + ")), data=dat)
    beta0 <- glht(model=fit1, linfct=paste("Intercept == 0"))
    beta1 <- glht(model=fit1, linfct=paste("Intercept + x*1 == 0"))
    beta2 <- glht(model=fit1, linfct=paste("Intercept + x*2 == 0"))
    varbeta0 <- glht(model=fit2, linfct=paste("Intercept == 0"))
    varbeta1 <- glht(model=fit2, linfct=paste("Intercept + x*1 + xsq*1 == 0"))
    varbeta2 <- glht(model=fit2, linfct=paste("Intercept + x*2 + xsq*4 == 0"))
    varbeta1_est <- glht(model=fit2, linfct=paste("x*1 + xsq*1 == 0"))
    varbeta2_est <- glht(model=fit2, linfct=paste("x*2 + xsq*4 == 0"))
    res <- cbind(
        tidy(beta0) %>% dplyr::select("estimate", "std.error") %>% dplyr::rename(estimate.0="estimate", std.error.0="std.error"),
        tidy(beta1) %>% dplyr::select("estimate", "std.error") %>% dplyr::rename(estimate.1="estimate", std.error.1="std.error"),
        tidy(beta2) %>% dplyr::select("estimate", "std.error") %>% dplyr::rename(estimate.2="estimate", std.error.2="std.error"),
        tidy(varbeta0) %>% dplyr::select("estimate", "std.error") %>% dplyr::rename(var.estimate.0="estimate", var.std.error.0="std.error"),
        tidy(varbeta1) %>% dplyr::select("estimate", "std.error") %>% dplyr::rename(var.estimate.1="estimate", var.std.error.1="std.error"),
        tidy(varbeta2) %>% dplyr::select("estimate", "std.error") %>% dplyr::rename(var.estimate.2="estimate", var.std.error.2="std.error"),
        tidy(varbeta1_est) %>% dplyr::select("estimate", "std.error") %>% dplyr::rename(var_effect.estimate.1="estimate", var_effect.std.error.1="std.error"),
        tidy(varbeta2_est) %>% dplyr::select("estimate", "std.error") %>% dplyr::rename(var_effect.estimate.2="estimate", var_effect.std.error.2="std.error"),
        n0=table(round(dat$x))[1],
        n1=table(round(dat$x))[2],
        n2=table(round(dat$x))[3],
        p=tidy(anova(fitnull, fit2))$p.value[2]
    )
    res$snp <- snp
    res$outcome <- outcome
    return(res)
}

get_lci <- function(b, se){
    b <- as.numeric(b)
    se <- as.numeric(se)
    lci <- b - (1.96 * se)
    return(lci)
}
get_uci <- function(b, se){
    b <- as.numeric(b)
    se <- as.numeric(se)
    uci <- b + (1.96 * se)
    return(uci)
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
d$id <- paste0("ukb-d-", str_split(d$trait, "\\.", simplify=T)[,2], "_irnt")
d$chr_pos <- paste0(d$chr, ":", d$pos)

# load finemapped snps
finemap <- fread("data/finemap.txt")
finemap.loci <- data.frame()
for (i in 1:nrow(finemap)){
    assoc <- ieugwasr::associations(variants=finemap$finemap.snp[i], id=finemap$id[i])
    finemap.loci <- rbind(finemap.loci, data.frame(
        finemap[i],
        chr=assoc$chr,
        pos=assoc$position,
        oa=assoc$ea,
        ea=assoc$nea,
        stringsAsFactors=F
    ), stringsAsFactors=F)
}
finemap.loci$key <- paste0("chr", finemap.loci$chr, "_", finemap.loci$pos, "_", finemap.loci$oa, "_", finemap.loci$ea)

# load dosage
snps <- d %>% dplyr::select(chr, pos, oa, ea)
snps <- rbind(snps, finemap.loci %>% dplyr::select(chr, pos, oa, ea), stringsAsFactors=F)
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
    
    # finemap variant analysis
    res_finemap <- bp(dat, d$key[i], d$trait[i], c("age_at_recruitment.21022.0.0", "sex.31.0.0", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", finemap.loci %>% dplyr::filter(trait==d$trait[i]) %>% dplyr::pull("key")))
    names(res_finemap) <- paste0(names(res_finemap), ".finemap")

    results <- rbind(results, cbind(res, res_log, res_finemap))
}

# write to table
write.table(results, file="data/vqtls.validate.txt", sep="\t", quote=F, row.names=F)

# write paper table
results$Trait <- NA
for (i in 1:nrow(results)){
    results$Trait[i] <- biomarkers_abr[biomarkers==results$outcome[i]]
}

results$var_effect.lci.1 <- get_lci(results$var_effect.estimate.1, results$var_effect.std.error.1)
results$var_effect.uci.1 <- get_uci(results$var_effect.estimate.1, results$var_effect.std.error.1)

results$var_effect.lci.1.log <- get_lci(results$var_effect.estimate.1.log, results$var_effect.std.error.1.log)
results$var_effect.uci.1.log <- get_uci(results$var_effect.estimate.1.log, results$var_effect.std.error.1.log)

results$var_effect.lci.1.finemap <- get_lci(results$var_effect.estimate.1.finemap, results$var_effect.std.error.1.finemap)
results$var_effect.uci.1.finemap <- get_uci(results$var_effect.estimate.1.finemap, results$var_effect.std.error.1.finemap)

results$var_effect.lci.2 <- get_lci(results$var_effect.estimate.2, results$var_effect.std.error.2)
results$var_effect.uci.2 <- get_uci(results$var_effect.estimate.2, results$var_effect.std.error.2)

results$var_effect.lci.2.log <- get_lci(results$var_effect.estimate.2.log, results$var_effect.std.error.2.log)
results$var_effect.uci.2.log <- get_uci(results$var_effect.estimate.2.log, results$var_effect.std.error.2.log)

results$var_effect.lci.2.finemap <- get_lci(results$var_effect.estimate.2.finemap, results$var_effect.std.error.2.finemap)
results$var_effect.uci.2.finemap <- get_uci(results$var_effect.estimate.2.finemap, results$var_effect.std.error.2.finemap)

tbl <- results %>% 
    dplyr::select(Trait, snp, var_effect.estimate.1, var_effect.lci.1, var_effect.uci.1, var_effect.estimate.2, var_effect.lci.2, var_effect.uci.2, p, var_effect.estimate.1.log, var_effect.lci.1.log, var_effect.uci.1.log, var_effect.estimate.2.log, var_effect.lci.2.log, var_effect.uci.2.log, p.log, var_effect.estimate.1.finemap, var_effect.lci.1.finemap, var_effect.uci.1.finemap, var_effect.estimate.2.finemap, var_effect.lci.2.finemap, var_effect.uci.2.finemap, p.finemap) %>% 
    dplyr::rename(beta.1="var_effect.estimate.1", lci.1="var_effect.lci.1", uci.1="var_effect.uci.1", beta.2="var_effect.estimate.2", lci.2="var_effect.lci.2", uci.2="var_effect.uci.2", beta.1.log="var_effect.estimate.1.log", lci.1.log="var_effect.lci.1.log", uci.1.log="var_effect.uci.1.log", beta.2.log="var_effect.estimate.2.log", lci.2.log="var_effect.lci.2.log", uci.2.log="var_effect.uci.2.log", beta.1.finemap="var_effect.estimate.1.finemap", lci.1.finemap="var_effect.lci.1.finemap", uci.1.finemap="var_effect.uci.1.finemap", beta.2.finemap="var_effect.estimate.2.finemap", lci.2.finemap="var_effect.lci.2.finemap", uci.2.finemap="var_effect.uci.2.finemap")

### add gene info ###

# read coloc data
coloc <- fread("data/coloc.txt")

# select high prob of shared causal variant
coloc <- coloc[coloc$PP.H4.abf > .8]

# append opengwas trait name to coloc data
ao <- available_outcomes()
ao <- ao %>% dplyr::select(c("trait", "id"))
coloc <- merge(coloc, ao, by.y="id", by.x="gene")

# append sun et al target info to coloc data
lookup <- fread("data/001_SOMALOGIC_GWAS_protein_info.csv")
lookup <- lookup %>% dplyr::select("Target", "TargetFullName")
coloc <- merge(coloc, lookup, by.x="trait.y", by.y="TargetFullName", all.x=T)

# append gene name to ENSG targets
bed <- fread("data/Homo_sapiens.GRCh37.82.bed")
bed <- bed %>% dplyr::select("V4", "V5")
coloc <- merge(coloc, bed, by.x="trait.y", by.y="V4", all.x=T)

# carry over eQTL gene names
coloc$Target <- apply(coloc,1,function(x) if( is.na(x[["Target"]]) ) x[["V5"]] else x[["Target"]]   )

# manually add in Folkersen et al targets
coloc[which(coloc$trait.y=="KIT ligand"),]$Target <- "KITLG"
coloc[which(coloc$trait.y=="chitinase 3 like 1"),]$Target <- "CHI3L1"
coloc[which(coloc$trait.y=="coagulation factor III, tissue factor"),]$Target <- "F3"
coloc[which(coloc$trait.y=="follistatin"),]$Target <- "FST"
coloc[which(coloc$trait.y=="galectin 3"),]$Target <- "LGALS3"
coloc[which(coloc$trait.y=="platelet and endothelial cell adhesion molecule 1"),]$Target <- "PECAM1"
coloc[which(coloc$trait.y=="selectin E"),]$Target <- "SELE"

# remove commas from gene name
coloc$Target <- gsub(",", "", coloc$Target)
coloc$Target <- gsub("  ", " ", coloc$Target)

# add neartest gene
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

s <- str_split(tbl$snp, "_", simplify=T)[,1:2]
s <- data.frame(chr=s[1], pos=s[2])
tbl <- cbind(tbl, s)
tbl$chr <- gsub("chr", "", tbl$chr)
tbl$key <- paste0(tbl$chr, ":", tbl$pos)

# merge
all <- merge(tbl, ng, "key")
all$key <- paste0(all$chr, ":", all$pos, "-", all$pos)
all$coloc <- NA

# add trait col
coloc <- coloc %>% filter(trait.x != "body_mass_index.21001.0.0")
coloc$Trait <- NA
for (i in 1:nrow(coloc)){
    coloc$Trait[i] <- biomarkers_abr[biomarkers==coloc$trait.x[i]]
}

# copy over
for (i in 1:nrow(all)){
    # snp region
    sregion <- GRanges(all$key[i])

    # gene region
    trait_coloc <- coloc %>% dplyr::filter(Trait == all$Trait[i])
    gregion <- GRanges(trait_coloc %>% dplyr::pull(region))

    # count overlap
    trait_coloc$counts <- suppressWarnings(countOverlaps(gregion, sregion))
    trait_coloc <- trait_coloc %>% dplyr::filter(counts >0 ) %>% unique()

    if (nrow(trait_coloc) == 0){
        next
    }
    
    # add gene with coloc
    trait_coloc <- trait_coloc %>% dplyr::select(Target, PP.H4.abf) %>% dplyr::group_by(Target) %>% dplyr::arrange(desc(PP.H4.abf)) %>% dplyr::filter(row_number()==1)
    trait_coloc$gene <- paste0(trait_coloc$Target, " (H4=", round(trait_coloc$PP.H4.abf,2), ")")
    trait_coloc <- trait_coloc %>% dplyr::arrange(desc(PP.H4.abf)) %>% head(n=3)
    f <- paste0(trait_coloc$gene, collapse="|")

    # set field in maste table
    all$coloc[i] <- f
}

all$key <- NULL
all$chr <- NULL
all$pos <- NULL
names(all)[names(all)=="gene"] <- "Nearest gene"
names(all)[names(all)=="coloc"] <- "Colocalised gene"
write.csv(all, file="data/supplementary_data_file_1.csv", row.names=F, quote=F)