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
d <- fread("data/gxg-qual.txt")
d$rsid.1 <- sapply(d$term, get_rsid)
d$rsid.2 <- sapply(d$mod_snp, get_rsid)
d$gene.1 <- sapply(d$term, get_gene)
d$gene.2 <- sapply(d$mod_snp, get_gene)
d$lci <- d$estimate - (1.96 * d$std.error)
d$uci <- d$estimate + (1.96 * d$std.error)

p1 <- ggplot(d, aes(x=mod, y=estimate, ymin=lci, ymax=uci, group=trait, shape=trait)) +
    geom_point(size = 1.5, position = position_dodge(width = 0.9)) +
    geom_errorbar(width=.05, position = position_dodge(width = 0.9)) +
    theme_classic() +
    facet_grid(Trait~., scales="free", space="free_y") + 
    scale_y_continuous(breaks = scales::pretty_breaks(5)) +
    geom_hline(yintercept = c(0), linetype = "dashed", color = "grey") +
    labs(shape = "Trait") +
    theme(
        panel.spacing.y = unit(0, "lines"),
        legend.box.background = element_rect(colour = "black"),
        legend.position = "none"
    ) +
    ylab("Trait (95% CI)") +
    xlab("SNP modifier (genotype copies)")

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