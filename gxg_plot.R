library("data.table")
library("ggplot2")
library("dplyr")
library("stringr")
library("ggpubr")
library("viridis")
library('forestplot')
library("RColorBrewer")
source("funs.R")
set.seed(123)

#read in gxg results
d <- fread("data/gxg.txt")
d$lci <- d$estimate - (d$std.error * 1.96)
d$uci <- d$estimate + (d$std.error * 1.96)

# read in WF
wf <- fread("data/gxg-wf.txt")
wf$lci <- wf$beta - (wf$se * 1.96)
wf$uci <- wf$beta + (wf$se * 1.96)

# filter genome-wide sig hits
d <- d %>% filter(p.value < 5e-8)

# merge
d <- merge(d, wf, "term")
d$t <- sapply(d$trait.x, function(x) biomarkers_abr[biomarkers==x]) %>% as.vector

# add rsid
lookup <- data.frame(
    snp=c(
        "chr11_116648917_C_G:chr8_19912370_A_G", 
        "chr12_121420260_G_A:chr19_45411941_C_T", 
        "chr22_44324730_T_C:chr4_88212722_A_G"
    ), 
    key=c(
        "APOA5 (rs964184G) x\nLPL (rs115849089G)",
        "HNF1A (rs7979473A) x\nAPOE (rs429358T)",
        "PNPLA3 (rs738408C) x\nHSD17B13 (rs13141441G)"
    )
)
d <- merge(d, lookup, by.x="term", by.y="snp")

# split out main & WF
main <- d %>% select(key, estimate, lci.x, uci.x, t) %>% rename(lci="lci.x", uci="uci.x")
main$analysis <- "Main"
wf <- d %>% select(key, beta, lci.y, uci.y, t) %>% rename(lci="lci.y", uci="uci.y", estimate="beta")
wf$analysis <- "Within-family"
d <- rbind(wf, main)
d$analysis <- factor(d$analysis, levels=c("Main", "Within-family"))

pdf("gxg.pdf", height=4)
ggplot(d, aes(x=key, y=estimate, ymin=lci, ymax=uci, group=analysis, color=t, shape=analysis)) +
    coord_flip() +
    scale_colour_brewer(palette = "Set1") +
    facet_grid(t~., scales="free", space="free_y") +
    geom_point(position = position_dodge(width = -0.25)) +
    geom_errorbar(width=.05, position = position_dodge(width = -0.25)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    scale_y_continuous(limits = c(-.2, .2), breaks=c(-.2, 0, .2)) +
    theme_classic() +
    labs(color = "Outcome", shape = "Model") +
    theme(
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        strip.text.y = element_blank()
    ) +
    ylab("Genotype * genotype (dosage) interaction effect estimate, SD (95% CI)")
dev.off()