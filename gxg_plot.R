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

# initialize plot
png("gxg-forest.png", width=480, height=480)

# make plot
forestplot(
    d$t,
    boxsize = 0.05,
    xticks = c(-.25, -.1, 0, .1, .25),
    legend = c("Main", "Within-family"),
    mean = d[,c("estimate", "beta")],
    lower = d[,c("lci.x","lci.y")],
    upper = d[,c("uci.x", "uci.y")],
    xlab=paste0("GxG effect (SD [95% CI])"),
    col=fpColors(lines=c("darkred", "royalblue"), box=c("darkred", "royalblue")),
    txt_gp = fpTxtGp(ticks=gpar(cex=1), xlab=gpar(cex=1))
)

# save plot
dev.copy(png, "gxg-forest.png")
dev.off()