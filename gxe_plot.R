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

get_dat <- function(file){
    # read in gxe results
    d <- fread(file)
    d$lci <- d$estimate - (d$std.error * 1.96)
    d$uci <- d$estimate + (d$std.error * 1.96)

    # filter genome-wide sig hits
    d <- d %>% filter(trait != "body_mass_index.21001.0.0")
    d <- d %>% filter(p.value < 5e-8)

    # merge
    d <- cbind(d, as.data.frame(str_split(d$term, ":", simplify=T), stringsAsFactors=F))
    d$y <- sapply(d$trait, function(x) return(biomarkers_abr[biomarkers==x]))
    d$u <- sapply(d$V1, get_trait_name)
    d <- d %>% filter(u != "Estimated Fat Yesterday") %>% filter(u != "Estimated Total Sugars Yesterday")

    # map SNP to rsid & gene
    lookup <- fread("all.vqtls.txt")
    lookup$key <- paste0(lookup$key, "-", lookup$Trait)
    d$key <- paste0(d$V2, "-", d$y)
    d <- merge(d, lookup, "key")
    d$key <- NULL
    d$gene <- str_split(d[["Nearest Gene"]], ",", simplify=T)[,1]
    d$f <- as.factor(paste0(d$gene, " (", d$RSID, d$EA, ")"))
    d$u <- factor(d$u)
    levels(d$u) <- list(Age="Age At Recruitment", Sex="Sex", BMI="Body Mass Index", PA="Summed Minutes Activity", Alcohol="Alcohol Intake Frequency", Smoking="Smoking Status")

    return(d)
}

get_plot <- function(d){
    # create row key
    key <- data.frame(Trait=sort(unique(d$Trait)), stringsAsFactors=F)
    key$key <- row(key) %% 2
    d <- merge(d, key, "Trait")
    d$key <- factor(d$key)

    # count number of traits
    n_traits <- length(unique(d$Trait))

    # Create a data frame with the faceting variables
    # and some dummy data (that will be overwritten)
    tp <- data.frame()
    for (tr in unique(d$Trait)){
        tp <- rbind(tp, data.frame(
            Trait=rep(tr, length(unique(d$u))),
            fill=which(tr == unique(d$Trait)) %% 2,
            u=unique(d$u)
        ))
    }
    tp$fill <- as.factor(tp$fill)

    p <- ggplot(d, aes(x=f, y=estimate, ymin=lci, ymax=uci)) +
        coord_flip() +
        facet_grid(Trait~u, scales="free", space="free_y") +
        geom_point(size = 1.5) +
        geom_errorbar(width=.05) +
        geom_hline(yintercept = c(-0.05, 0, 0.05), linetype = "dashed", color = "grey") +
        theme_classic() +
        scale_y_continuous(limits = c(-.1, .1), breaks=c(-.1, 0, .1)) +
        geom_rect(inherit.aes = F, data = tp, aes(fill = fill), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.15) +
        scale_fill_manual(values=brewer.pal(2,"Paired")) +
        theme(
            axis.title.y = element_blank(),
            strip.background = element_blank(),
            strip.text.y = element_text(angle = 0),
            legend.position = "none",
            panel.spacing.y = unit(0, "lines")
        ) +
        ylab("Genotype (dosage) * modifier (SD) interaction effect estimate, SD (95% CI)")
    return(p)
}

# load gxe effects
d <- get_dat("data/gxe.txt")
d <- get_dat("data/gxe-log.txt")

# save plot
pdf("gxe.pdf", height=11, width=11)
print(get_plot(d))
dev.off()