library("ieugwasr")
library("dplyr")
library("multcomp")
library("broom")
library('robustbase')
library("data.table")
library("ggplot2")
library("dplyr")
library("stringr")
library("ggpubr")
library("viridis")
library('forestplot')
library("viridis")
library("RColorBrewer")
library("grid")
source("funs.R")
set.seed(1234)
options(ieugwasr_api="http://64.227.44.193:8006/")

# load phenotypes
load("data/pheno.RData")
linker <- get_filtered_linker(drop_standard_excl=TRUE, drop_non_white_british=TRUE, drop_related=TRUE, application="15825")
covariates <- get_covariates()
pc <- get_genetic_principal_components()
dat <- merge(linker, covariates, "appieu")
dat <- merge(dat, pheno, by.x="app15825", by.y="eid")
dat <- merge(dat, pc, "appieu")

# split env into binary groups
dat$smoking_status.20116.0.0_b <- dat$smoking_status.20116.0.0 > 0
dat$summed_minutes_activity.22034.0.0_b <- dat$summed_minutes_activity.22034.0.0 < median(dat$summed_minutes_activity.22034.0.0, na.rm=T)
dat$alcohol_intake_frequency.1558.0.0_b <- dat$alcohol_intake_frequency.1558.0.0 < 4
dat$estimated_fat_yesterday.100004.0.0_b <- dat$estimated_fat_yesterday.100004.0.0 > median(dat$estimated_fat_yesterday.100004.0.0, na.rm=T)
dat$estimated_total_sugars_yesterday.100008.0.0_b <- dat$estimated_total_sugars_yesterday.100008.0.0 > median(dat$estimated_total_sugars_yesterday.100008.0.0, na.rm=T)
dat$age_at_recruitment.21022.0.0_b <- dat$age_at_recruitment.21022.0.0 > median(dat$age_at_recruitment.21022.0.0, na.rm=T)
dat$body_mass_index.21001.0.0_b <- dat$body_mass_index.21001.0.0 > median(dat$body_mass_index.21001.0.0, na.rm=T)
dat$sex.31.0.0_b <- dat$sex.31.0.0 > 0

# SD outcome
dat$ldl_direct.30780.0.0 <- dat$ldl_direct.30780.0.0 / sd(dat$ldl_direct.30780.0.0, na.rm=T)
dat$urate.30880.0.0 <- dat$urate.30880.0.0 / sd(dat$urate.30880.0.0, na.rm=T)

# get eQTL instruments
lpl <- fread("data/lpl.txt")
lpl <- lpl %>% dplyr::filter(chr == "8" & position > 19759228-500000 & position < 19824769+500000)
lpl$gene <- "LPL"

slc2a9 <- fread("data/slc2a9.txt")
slc2a9 <- slc2a9 %>% dplyr::filter(chr == "4" & position > 9772777-500000 & position < 10056560+500000)
slc2a9$gene <- "SLC2A9"

# combine
snps <- rbind(lpl, slc2a9)
snps$V1 <- NULL

results <- data.frame()
for (i in 1:nrow(snps)){
    # load dosage
    dosage <- extract_variant_from_bgen(as.character(snps$chr[i]), as.double(snps$position[i]), snps$nea[i], snps$ea[i])
    names(dosage)[1] <- snps$rsid[i]
    dat <- merge(dat, dosage, "appieu")

    # test for SNP-outcome effect stratified by modifier
    if (snps$gene[i] == "LPL"){
        o <- "ldl_direct.30780.0.0"
        k <- "body_mass_index.21001.0.0_b" 
    } else if (snps$gene[i] == "SLC2A9"){
        o <- "urate.30880.0.0"
        k <- "sex.31.0.0_b"
    }

    f <- as.formula(paste0(o, " ~ ", snps$rsid[i], " + age_at_recruitment.21022.0.0 + sex.31.0.0 + PC1 + PC2 + PC3 + PC4 + PC6 + PC7 + PC8 + PC9 + PC10"))
    t <- tidy(lmrob(f, data=dat %>% dplyr::filter(!!sym(k) == T)))
    res1 <- t[2,]
    t <- tidy(lmrob(f, data=dat %>% dplyr::filter(!!sym(k) == F)))
    res2 <- t[2,]
    names(res1) <- paste0(names(res1), ".T")
    names(res2) <- paste0(names(res2), ".F")
    res <- cbind(res1, res2)
    res$u <- k
    res$y <- o
    results <- rbind(results, res)
}

# plot
# create row key
d <- results
d$Trait <- d$y
d$lci.T <- d$estimate.T - (d$std.error.T * 1.96)
d$uci.T <- d$estimate.T + (d$std.error.T * 1.96)
d$lci.F <- d$estimate.F - (d$std.error.F * 1.96)
d$uci.F <- d$estimate.F + (d$std.error.F * 1.96)
key <- data.frame(Trait=sort(unique(d$Trait)), stringsAsFactors=F)
key$key <- row(key) %% 2
d <- merge(d, key, "Trait")
d$key <- factor(d$key)

# count number of traits
n_traits <- length(unique(d$Trait))

# Create a data frame with the faceting variables
# and some dummy data (that will be overwritten)
tp <- data.frame()
for (tr in unique(d$Trait)){
    tp <- rbind(tp, data.frame(
        Trait=rep(tr, length(unique(d$u))),
        fill=which(tr == unique(d$Trait)) %% 2,
        u=unique(d$u)
    ))
}
tp$fill <- as.factor(tp$fill)

# wide to long
e <- d
f <- d

e$estimate.F <- NULL
e$std.error.F <- NULL
e$statistic.F <- NULL
e$p.value.F <- NULL
e$lci.F <- NULL
e$uci.F <- NULL
e$subgroup <- T
e <- e %>% dplyr::rename(estimate="estimate.T", std.error="std.error.T", statistic="statistic.T", p.value="p.value.T", lci="lci.T", uci="uci.T")

f$estimate.T <- NULL
f$std.error.T <- NULL
f$statistic.T <- NULL
f$p.value.T <- NULL
f$lci.T <- NULL
f$uci.T <- NULL
f$subgroup <- F
f <- f %>% dplyr::rename(estimate="estimate.F", std.error="std.error.F", statistic="statistic.F", p.value="p.value.F", lci="lci.F", uci="uci.F")

d <- rbind(e, f)

# create plot
d$subgroup <- as.character(d$subgroup)
d[which(d$subgroup == "FALSE"),]$subgroup <- "Q1"
d[which(d$subgroup == "TRUE"),]$subgroup <- "Q2"
p <- ggplot(d, aes(x=term.T, y=estimate, ymin=lci, ymax=uci, color=subgroup)) +
    coord_flip() +
    facet_grid(Trait~u, scales="free", space="free_y") +
    geom_point(size = 1.5, position = position_dodge(width = 0.9)) +
    geom_errorbar(width=.05, position = position_dodge(width = 0.9)) +
    theme_classic() +
    geom_rect(inherit.aes = F, show.legend = FALSE, data = tp, aes(fill = fill), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.15) +
    scale_y_continuous(breaks = scales::pretty_breaks(3)) +
    geom_hline(yintercept = c(0), linetype = "dashed", color = "grey") +
    scale_fill_manual(values=brewer.pal(2,"Paired")[1:2]) +
    scale_color_manual(values=brewer.pal(4,"Paired")[3:4], name="Subgroup") +
    theme(
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        strip.text.y = element_text(angle = 0),
        legend.position = "bottom",
        panel.spacing.y = unit(0, "lines"),
        legend.box.background = element_rect(colour = "black")
    ) +
    ylab("Per allele effect estimate stratified by modifier, SD (95% CI)")