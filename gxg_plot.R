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
d$t <- get_trait_name(d$trait.x)

# add rsid
lookup <- data.frame(
    snp=c(
        "chr11_116648917_C_G:chr8_19912370_A_G", 
        "chr12_121420260_G_A:chr19_45411941_C_T", 
        "chr22_44324730_T_C:chr4_88212722_A_G"
    ), 
    rsid=c(
        "rs964184-G x rs115849089-G",
        "rs7979473-A x rs429358-T",
        "rs738408-C x rs13141441-G"
    ),
    gene=c(
        "APOA5 x LPL",
        "HNF1A x APOE",
        "PNPLA3 x HSD17B13"
    )
)
d <- merge(d, lookup, by.x="term", by.y="snp")

# plot
pdf("gxg-forest.pdf")
forestplot(
    paste0(d$t, "\n", d$gene, "\n", d$rsid),
    boxsize = 0.05,
    xticks = c(-.25, -.1, 0, .1, .25),
    legend = c("Main", "Within-family"),
    mean = d[,c("estimate", "beta")],
    lower = d[,c("lci.x","lci.y")],
    upper = d[,c("uci.x", "uci.y")],
    new_page=F,
    xlab=paste0("GxG effect (SD [95% CI])"),
    col=fpColors(lines=c("darkred", "royalblue"), box=c("darkred", "royalblue")),
    txt_gp = fpTxtGp(ticks=gpar(cex=1), xlab=gpar(cex=1))
)
dev.off()