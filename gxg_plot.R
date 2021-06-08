library("data.table")
library("ggplot2")
library("dplyr")
library("stringr")
library("ggpubr")
library("viridis")
set.seed(123)

#read in gxg results
d <- fread("data/gxg.txt")
d$lci <- d$estimate - (d$std.error * 1.96)
d$uci <- d$estimate + (d$std.error * 1.96)

# drop oestradiol.30800.0.0 which excess t1e
d <- d %>% filter(trait!="oestradiol.30800.0.0")

p <- ggplot(d, aes(x=x, y=estimate, ymin=lci, ymax=uci, group=trait, color=trait)) +
    geom_point(position=position_dodge(width=0.75)) +
    geom_errorbar(width=.05, position=position_dodge(width=0.75)) +
    theme_classic() + 
    ggtitle("GxG effects") +
    xlab("Gene") + 
    ylab(paste0("Interaction effect (SD, 95% CI)")) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    labs(col="Trait") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

pdf("gxg.pdf")
print(p)
dev.off()