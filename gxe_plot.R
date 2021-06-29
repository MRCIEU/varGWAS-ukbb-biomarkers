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
set.seed(123)

# read in gxe results
d <- fread("data/gxe.txt")
d$lci <- d$estimate - (d$std.error * 1.96)
d$uci <- d$estimate + (d$std.error * 1.96)

# filter genome-wide sig hits
d <- d %>% filter(trait != "body_mass_index.21001.0.0")
d <- d %>% filter(p.value < 5e-8)

# merge
d <- cbind(d, as.data.frame(str_split(d$term, ":", simplify=T), stringsAsFactors=F))
d$y <- sapply(d$trait, function(x) return(biomarkers_abr[biomarkers==x]))
d$u <- sapply(d$V1, get_trait_name)
d <- d %>% filter(u != "Estimated Fat Yesterday") %>% filter(u != "Estimated Total Sugars Yesterday") %>% filter(u != "Summed Minutes Activity")

# map SNP to rsid & gene
lookup <- fread("all.vqtls.txt")
lookup$key <- paste0(lookup$key, "-", lookup$Trait)
d$key <- paste0(d$V2, "-", d$y)
d <- merge(d, lookup, "key")
d$key <- NULL
d$f <- as.factor(d$RSID)
d$u <- factor(d$u)
levels(d$u) <- list(Age="Age At Recruitment", Alcohol="Alcohol Intake Frequency", BMI="Body Mass Index", Sex="Sex", Smoking="Smoking Status", PA="Summed Minutes Activity")

# select fields
main <- d %>% select(Trait, f, estimate, lci, uci, u, p.value)
main$analysis <- "Main"
main$logP <- -log10(main$p.value)

ggplot(d, aes(x=f, y=estimate, ymin=lci, ymax=uci)) +
    coord_flip() +
    facet_grid(Trait~u, scales="free", space="free_y") +
    geom_point() +
    geom_errorbar(width=.05) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    theme_classic() +
    scale_y_continuous(breaks = scales::breaks_pretty(n=3)) +
    labs(color = expression(paste("-log"[10],"(", plain(P),")"))) +
    theme(
        strip.text.y = element_text(angle=0),
        axis.title.y = element_blank(),
        strip.background.y = element_blank()
    ) +
    ylab("Genotype (dosage) * modifier (SD) interaction effect estimate, SD (95% CI)")