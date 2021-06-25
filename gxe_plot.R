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
d <- fread("data/gxe.txt")
d$lci <- d$estimate - (d$std.error * 1.96)
d$uci <- d$estimate + (d$std.error * 1.96)

# read in WF
wf <- fread("data/gxe-wf.txt")
wf$lci <- wf$beta - (wf$se * 1.96)
wf$uci <- wf$beta + (wf$se * 1.96)
wf$term_trait <- paste0(wf$term, "-", wf$trait)

# filter genome-wide sig hits
d <- d %>% filter(trait != "body_mass_index.21001.0.0")
d <- d %>% filter(p.value < 5e-8)
d$term_trait <- paste0(d$term, "-", d$trait)

# merge
d <- merge(d, wf, "term_trait")
d$t <- get_trait_name(d$trait.x)

# plot
pdf("gxe-forest.pdf")
forestplot(
    d$term_trait,
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