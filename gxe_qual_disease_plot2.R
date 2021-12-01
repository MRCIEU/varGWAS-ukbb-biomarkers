load("data/pheno.RData")
library('optparse')
library('data.table')
library('dplyr')
library('broom')
library('ieugwasr')
library("ggpubr")
library("lmtest")
library("ggplot2")
library("sandwich")
library("susieR")
source("funs.R")
set.seed(1234)
options(ieugwasr_api="http://64.227.44.193:8006/")

finemap <- function(dat){
    # format data
    dat <- dat %>% filter(grepl("^rs", term))
    dat <- cbind(dat, str_split(dat$term, ":", simplify=T)) %>% rename(SNP="V1",modifier="V2")
    dat <- cbind(dat, str_split(dat$SNP, "_", simplify=T)) %>% rename(rsid="V1", oa="V2", ea="V3")
    dat$term <- NULL
    dat$statistic <- NULL
    dat$zscore <- dat$estimate / dat$std.error
    
    ld <- ieugwasr::ld_matrix (
        dat$rsid,
        with_alleles = T,
        pop = "EUR",
        bfile = "/user/home/ml18692/varGWAS-ukbb-biomarkers/data/EUR",
        plink_bin = "/user/home/ml18692/varGWAS-ukbb-biomarkers/data/plink_Linux"
    )

    # select SNPs with LD data and orient effect size
    dat2 <- dat %>% filter(SNP %in% row.names(ld))
    dat$SNP <- paste0(dat$rsid, "_", dat$ea, "_", dat$oa)
    dat$estimate <- dat$estimate *-1
    dat$zscore <- dat$estimate / dat$std.error
    dat3 <- dat %>% filter(SNP %in% row.names(ld))
    dat <- rbind(dat2, dat3)
    dat$rsid <- NULL
    dat$oa <- NULL
    dat$ea <- NULL

    # order by LD matrix
    dat4 <- data.frame()
    for (SNP in row.names(ld)){
        dat4 <- rbind(dat4, dat %>% filter(SNP == !!SNP))
    }
    
    fitted_rss <- tryCatch(
        expr = {
            susieR::susie_rss(
                dat4$zscore,
                ld,
                L=10,
                estimate_prior_variance=TRUE,
                max_iter=500
            )
        },
        error = function(e){ 
            NULL
        }
    )

    if (is.null(fitted_rss)){
        return(NULL)
    }

    # collect fine mapped snps
    cs <- summary(fitted_rss)$cs

    if (is.null(cs)){
        return(NULL)
    }

    snps <- data.frame()
    for (j in 1:nrow(cs)){
        for (k in stringr::str_split(cs$variable[j], ",", simplify=T) %>% as.numeric){
            rsid=str_split(dat4$SNP[k], "_", simplify=T)[1]
            res <- cs[j,]
            res <- cbind(res, data.frame(
                snp=dat4$SNP[k],
                pip=fitted_rss$pip[k]
            ))
            res <- cbind(res, associations(rsid, "ukb-d-30640_irnt") %>% select(chr, position, nea, ea))
            snps <- rbind(snps, res)
        }
    }

    return(snps)
}

get_lm_estimate <- function(f, k, dat, trait, g1, g2){
    mod <- lm(f, data=dat %>% dplyr::filter(!!sym(k) == T))
    t <- coeftest(mod, vcov = vcovHC(mod, type = "HC0")) %>% tidy
    res1 <- t[2,]
    mod <- lm(f, data=dat %>% dplyr::filter(!!sym(k) == F))
    t <- coeftest(mod, vcov = vcovHC(mod, type = "HC0")) %>% tidy
    res2 <- t[2,]
    mod <- lm(f, data=dat)
    t <- coeftest(mod, vcov = vcovHC(mod, type = "HC0")) %>% tidy
    res3 <- t[2,]
    res1$group <- g1
    res2$group <- g2
    res3$group <- "All"
    res <- rbind(res1, res2, res3)
    res$trait <- trait
    res$binary <- FALSE
    res$k <- k
    return(res)
}

get_glm_estimate <- function(f, k, dat, trait, g1, g2){
    mod <- glm(f, family="binomial", data=dat %>% dplyr::filter(!!sym(k) == T))
    t <- tidy(mod)
    res1 <- t[2,]
    mod <- glm(f, family="binomial", data=dat %>% dplyr::filter(!!sym(k) == F))
    t <- tidy(mod)
    res2 <- t[2,]
    mod <- glm(f, family="binomial", data=dat)
    t <- tidy(mod)
    res3 <- t[2,]
    res1$group <- g1
    res2$group <- g2
    res3$group <- "All"
    res <- rbind(res1, res2, res3)
    res$trait <- trait
    res$binary <- TRUE
    res$k <- k
    return(res)
}

est <- function(f_dat, vqtl, dat, trait_c, trait_b, group, k1, k2){
    results <- data.frame()
    chrsnp <- paste0("chr", vqtl$chr, "_", vqtl$position, "_", vqtl$nea, "_", vqtl$ea)
    lookup <- data.frame(rsid=str_split(f_dat$snp, "_", simplify=T)[,1], key=paste0("chr", f_dat$chr, "_", f_dat$position, "_", f_dat$nea, "_", f_dat$ea), stringsAsFactors=F)
    lookup <- rbind(lookup, data.frame(rsid=vqtl$rsid, key=chrsnp, stringsAsFactors=F))
    lookup <- unique(lookup)

    # load finemap dosage
    for (i in 1:nrow(f_dat)){
        dosage <- extract_variant_from_bgen(f_dat$chr[i], as.double(f_dat$position[i]), f_dat$nea[i], f_dat$ea[i])
        dat <- merge(dat, dosage, "appieu")
    }
    
    # load vQTL dosage
    if (!chrsnp %in% names(dat)){
        dosage <- extract_variant_from_bgen(vqtl$chr, as.double(vqtl$position), vqtl$nea, vqtl$ea)
        dat <- merge(dat, dosage, "appieu")
    }

    for (snp in unique(c(chrsnp, paste0("chr", f_dat$chr, "_", f_dat$position, "_", f_dat$nea, "_", f_dat$ea)))){
        f <- as.formula(paste0(trait_c, " ~ ", snp, " + age_at_recruitment.21022.0.0 + sex.31.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))
        results <- rbind(results, get_lm_estimate(f, group, dat, trait_c, k1, k2))
        f <- as.formula(paste0(trait_b, " ~ ", snp, " + age_at_recruitment.21022.0.0 + sex.31.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))
        results <- rbind(results, get_glm_estimate(f, group, dat, trait_b, k1, k2))
    }

    results$lci <- results$estimate - (1.96 * results$std.error)
    results$uci <- results$estimate + (1.96 * results$std.error)
    results$estimate[results$binary] <- exp(results$estimate[results$binary])
    results$lci[results$binary] <- exp(results$lci[results$binary])
    results$uci[results$binary] <- exp(results$uci[results$binary])

    results <- merge(results, lookup, by.x="term", by.y="key")
    #results$rsid <- factor(results$rsid, levels=unique(c(vqtl$rsid, f_dat$rsid)))

    f_dat$rsid <- str_split(f_dat$snp, "_", simplify=T)[,1]
    results <- merge(results, f_dat %>% select(cs, rsid), "rsid", all.x=T)
    results$cs[which(results$rsid == vqtl$rsid)] <- 0
    results <- results %>% arrange(cs)
    results$cs <- paste0("CS ", results$cs)
    results$cs <- gsub("CS 0", "vQTL", results$cs)
    results$cs <- factor(results$cs, levels=c("vQTL", "CS 1", "CS 2", "CS 3", "CS 4", "CS 5", "CS 6", "CS 7", "CS 8", "CS 9"))
        
    return(results)
}

plot <- function(dat_est, xlab, ylab){
    p <- ggplot(dat_est, aes(x=group, y=estimate, ymin=lci, ymax=uci, color=rsid)) +
        geom_point(size = 1.5, position = position_dodge(width = 0.9)) +
        geom_errorbar(width=.05, position = position_dodge(width = 0.9)) +
        theme_classic() +
        facet_grid(~ cs) +
        scale_y_continuous(breaks = scales::pretty_breaks(5)) +
        geom_hline(yintercept = c(0), linetype = "dashed", color = "grey") +
        theme(
            panel.spacing.y = unit(0, "lines"),
            legend.box.background = element_rect(colour = "black"),
            legend.position = "bottom",
            legend.title = element_blank()
        ) +
        xlab(xlab) +
        ylab(ylab)
    return(p)
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

# split env into binary groups
dat$body_mass_index.21001.0.0_b <- dat$body_mass_index.21001.0.0 > median(dat$body_mass_index.21001.0.0, na.rm=T)
dat$sex.31.0.0_b <- dat$sex.31.0.0 > 0

# SD outcome
dat$hdl_cholesterol.30760.0.0 <- dat$hdl_cholesterol.30760.0.0 / sd(dat$hdl_cholesterol.30760.0.0, na.rm=T)
dat$triglycerides.30870.0.0 <- dat$triglycerides.30870.0.0 / sd(dat$triglycerides.30870.0.0, na.rm=T)
dat$urate.30880.0.0 <- dat$urate.30880.0.0 / sd(dat$urate.30880.0.0, na.rm=T)

# load GxE interval around vQTL
rs75627662 <- fread("data/hdl_cholesterol.30760.0.0.x.sex.31.0.0.rs75627662.txt") %>% filter(grepl(":", term))
rs738409 <- fread("data/triglycerides.30870.0.0.x.body_mass_index.21001.0.0.rs738409.txt") %>% filter(grepl(":", term))
rs4530622 <- fread("data/urate.30880.0.0.x.sex.31.0.0.rs4530622.txt2") %>% filter(grepl(":", term))

# finemap GxE effects
rs75627662_f <- finemap(rs75627662)
rs738409_f <- finemap(rs738409)
rs4530622_f <- finemap(rs4530622)

# estimate per SNP effect on outcome by subgroup of vQTL and finemapped variants
rs75627662_est <- est(rs75627662_f, data.frame(chr="19", position=45413576, nea="C", ea="T", rsid="rs75627662", stringsAsFactors=F), dat, "hdl_cholesterol.30760.0.0", "vascular_problems.6150", "sex.31.0.0_b", "Male", "Female")
rs738409_est <- est(rs738409_f, data.frame(chr="22", position=44324727, nea="C", ea="G", rsid="rs738409", stringsAsFactors=F), dat, "triglycerides.30870.0.0", "vascular_problems.6150", "body_mass_index.21001.0.0_b", "BMI > 26.7 kg/m2", "BMI < 26.7 kg/m2")
rs4530622_est <- est(rs4530622_f, data.frame(chr="4", position=10402838, nea="T", ea="C", rsid="rs4530622", stringsAsFactors=F), dat, "urate.30880.0.0", "gout", "sex.31.0.0_b", "Male", "Female")

rs4530622_est$cs[which(rs4530622_est$cs == "CS 6")] <- "CS 5"
rs4530622_est$cs[which(rs4530622_est$cs == "CS 7")] <- "CS 6"
rs4530622_est$cs[which(rs4530622_est$cs == "CS 8")] <- "CS 7"
rs4530622_est$cs[which(rs4530622_est$cs == "CS 9")] <- "CS 8"
rs4530622_est$cs <- factor(as.character(rs4530622_est$cs), levels=c("vQTL", "CS 1", "CS 2", "CS 3", "CS 4", "CS 5", "CS 6", "CS 7", "CS 8"))

p1 <- plot(rs75627662_est %>% filter(trait == "hdl_cholesterol.30760.0.0"), xlab = "Sex", ylab = "HDL (SD, 95% CI)")
p2 <- plot(rs738409_est %>% filter(trait == "triglycerides.30870.0.0"), xlab = "BMI", ylab = "TG (SD, 95% CI)")
p3 <- plot(rs4530622_est %>% filter(trait == "urate.30880.0.0"), xlab = "Sex", ylab = "Urate (SD, 95% CI)")

pa <- ggarrange(p1, p2, labels = c("A", "B"), ncol = 2, nrow = 1)
pb <- ggarrange(pa, p3, labels = c(NA, "C"), ncol = 1, nrow = 2)

pdf("gxe-qual2.pdf", height=7*2, width=16)
print(pb)
dev.off()