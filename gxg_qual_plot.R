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

# read in gxg qual
d <- fread("data/alanine_aminotransferase.30620.0.0.gxg-qual.txt")
d$lci <- d$estimate - (1.96 * d$std.error)
d$uci <- d$estimate + (1.96 * d$std.error)
d$mod_gt <- NA
d$mod_gt[d$mod == 0] <- "TT"
d$mod_gt[d$mod == 1] <- "TC"
d$mod_gt[d$mod == 2] <- "CC"

p <- ggplot(d, aes(x=mod_gt, y=estimate, ymin=lci, ymax=uci)) +
    geom_point(size = 1.5, position = position_dodge(width = 0.9)) +
    geom_errorbar(width=.05, position = position_dodge(width = 0.9)) +
    theme_classic() +
    scale_y_continuous(breaks = scales::pretty_breaks(5)) +
    geom_hline(yintercept = c(0), linetype = "dashed", color = "grey") +
    theme(
        panel.spacing.y = unit(0, "lines"),
        legend.box.background = element_rect(colour = "black")
    ) +
    ylab("Per allele effect estimate, (95% CI)") +
    xlab("PNPLA3 (rs738408)")

pdf("gxg-qual.pdf")
print(p)
dev.off()