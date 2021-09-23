library("data.table")
library("ggplot2")
library("dplyr")
library("stringr")
library("ggpubr")
library("viridis")
library('forestplot')
library("viridis")
library("ggpubr")
library("RColorBrewer")
library("grid")
source("funs.R")
set.seed(123)

# read in gxg qual
d <- fread("data/alanine_aminotransferase.30620.0.0.gxg-qual.txt")
d$lci <- d$estimate - (1.96 * d$std.error)
d$uci <- d$estimate + (1.96 * d$std.error)
d$mod_gt <- NA
d$mod_gt[d$mod == 0] <- "GG"
d$mod_gt[d$mod == 1] <- "GA"
d$mod_gt[d$mod == 2] <- "AA"
d$mod_gt <- factor(d$mod_gt, levels=c("GG", "GA", "AA"))

d1 <- d %>% filter(trait=="alanine_aminotransferase.30620.0.0")
d1$trait <- gsub("alanine_aminotransferase.30620.0.0", "Alanine aminotransferase", d1$trait)
d2 <- d %>% filter(trait!="alanine_aminotransferase.30620.0.0")
d2 <- d %>% filter(trait=="liver_disease")
d2$estimate <- exp(d2$estimate)
d2$lci <- exp(d2$lci)
d2$uci <- exp(d2$uci)
d2$trait <- gsub("alcoholic_liver_disease", "Alcoholic liver disease", d2$trait)
d2$trait <- gsub("fibrosis_liver_disease", "Fibrotic liver disease", d2$trait)
d2$trait <- gsub("fatty_liver_disease", "Fatty liver disease", d2$trait)
d2$trait <- gsub("liver_disease", "Liver disease", d2$trait)

p1 <- ggplot(d1, aes(x=mod_gt, y=estimate, ymin=lci, ymax=uci, group=trait, shape=trait)) +
    geom_point(size = 1.5, position = position_dodge(width = 0.9)) +
    geom_errorbar(width=.05, position = position_dodge(width = 0.9)) +
    theme_classic() +
    scale_y_continuous(breaks = scales::pretty_breaks(5)) +
    geom_hline(yintercept = c(0), linetype = "dashed", color = "grey") +
    labs(shape = "Trait") +
    theme(
        panel.spacing.y = unit(0, "lines"),
        legend.box.background = element_rect(colour = "black"),
        legend.position = "none"
    ) +
    ylab("ALT (95% CI)") +
    xlab("HSD17B13 (rs13141441)")

p2 <- ggplot(d2, aes(x=mod_gt, y=estimate, ymin=lci, ymax=uci, group=trait, shape=trait)) +
    geom_point(size = 1.5, position = position_dodge(width = 0.9)) +
    geom_errorbar(width=.05, position = position_dodge(width = 0.9)) +
    theme_classic() +
    scale_y_continuous(breaks = scales::pretty_breaks(5)) +
    geom_hline(yintercept = c(1), linetype = "dashed", color = "grey") +
    labs(shape = "Trait") +
    theme(
        panel.spacing.y = unit(0, "lines"),
        legend.box.background = element_rect(colour = "black"),
        legend.position = "none"
    ) +
    ylab("Liver disease (OR, 95% CI)") +
    xlab("HSD17B13 (rs13141441)")

p <- ggarrange(p1, p2, labels = c("A", "B"), ncol = 2, nrow = 1)

pdf("gxg-qual.pdf", height=7*(2/3), width=14)
print(p)
dev.off()