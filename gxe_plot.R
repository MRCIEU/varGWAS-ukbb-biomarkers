library("data.table")
library("ggplot2")
library("dplyr")
library("stringr")
library("ggpubr")
library("viridis")
library('forestplot')
library("RColorBrewer")
library("grid")
source("funs.R")
set.seed(123)

#read in gxg results
d <- fread("data/gxe.txt")
d$lci <- d$estimate - (d$std.error * 1.96)
d$uci <- d$estimate + (d$std.error * 1.96)

# read in WF
wf <- fread("data/gxe-wf.txt")
wf$lci <- wf$beta - (wf$se * 1.96)
wf$uci <- wf$beta + (wf$se * 1.96)
wf$key <- paste0(wf$term, "-", wf$trait)
wf$trait <- NULL
wf$term <- NULL

# filter genome-wide sig hits
d <- d %>% filter(trait != "body_mass_index.21001.0.0")
d <- d %>% filter(p.value < 5e-8)
d$key <- paste0(d$term, "-", d$trait)

# merge
d <- merge(d, wf, "key")
d$key <- NULL
d <- cbind(d, as.data.frame(str_split(d$term, ":", simplify=T), stringsAsFactors=F))
d$y <- sapply(d$trait, function(x) return(biomarkers_abr[biomarkers==x]))
d$u <- sapply(d$V1, get_trait_name)

# map SNP to rsid & gene
lookup <- fread("all.vqtls.txt")
lookup$key <- paste0(lookup$key, "-", lookup$Trait)
d$key <- paste0(d$V2, "-", d$y)
d <- merge(d, lookup, "key")
d$key <- NULL
d$f <- as.factor(d$RSID)

# plot
d$u <- factor(d$u)
levels(d$u) <- list(Age="Age At Recruitment", Smoking="Smoking Status", Sex="Sex", BMI="Body Mass Index", Alcohol="Alcohol Intake Frequency")
ggplot(data=d, aes(x=f, y=estimate, ymin=lci.x, ymax=uci.x, group=Trait, color=Trait)) +
    coord_flip() +
    facet_grid(u~., scales = "free_y", space="free_y")+
    geom_point() +
    geom_errorbar(width=.05) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    theme_classic() +
    xlab("SNP") + 
    ylab("Estimate, SD (95% CI)") +
    theme(
        strip.background = element_blank(),
        strip.text.y = element_text(angle=0),
        panel.spacing = unit(2, "lines")
    )

ggplot(data=d, aes(x=f, y=estimate, ymin=lci.x, ymax=uci.x, group=u, color=u)) +
    coord_flip() +
    facet_grid(Trait~., scales = "free_y", space="free_y")+
    geom_point() +
    geom_errorbar(width=.05) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    theme_classic() +
    xlab("SNP") + 
    ylab("Estimate, SD (95% CI)") +
    theme(
        strip.text.y = element_text(angle=0),
        strip.background = element_blank()
    )