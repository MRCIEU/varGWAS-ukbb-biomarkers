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
    d$term.F <- NULL
    names(d)[1] <- "term"
    d$lci.T <- d$estimate.T - (d$std.error.T * 1.96)
    d$uci.T <- d$estimate.T + (d$std.error.T * 1.96)
    d$lci.F <- d$estimate.F - (d$std.error.F * 1.96)
    d$uci.F <- d$estimate.F + (d$std.error.F * 1.96)
    d$key <- paste0(d$term, ":", d$u, ":", d$y)

    # drop BMI
    d <- d %>% filter(y != "body_mass_index.21001.0.0")

    # filter SNPs to show
    int <- fread("data/gxg.txt")
    int$key <- paste0(int$term, ":", int$trait)
    int <- int %>% dplyr::filter(p.value < 5e-8)
    int_log <- fread("data/gxg-log.txt")
    int_log$key <- paste0(int_log$term, ":", int_log$trait)
    int_log <- int_log %>% dplyr::filter(p.value < 5e-5)
    int <- int %>% dplyr::filter(key %in% int_log$key)
    d$key_mb <- gsub("_bin", "", d$key)
    d <- d %>% dplyr::filter(key_mb %in% int$key)
    d$y <- sapply(d$y, function(x) return(biomarkers_abr[biomarkers==x]))
    d$key <- NULL
    d$key_mb <- NULL
    d$Trait <- d$y

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

    # wide to long
    e <- d
    f <- d

    e$estimate.F <- NULL
    e$std.error.F <- NULL
    e$statistic.F <- NULL
    e$p.value.F <- NULL
    e$lci.F <- NULL
    e$uci.F <- NULL
    e$subgroup <- T
    e <- e %>% dplyr::rename(estimate="estimate.T", std.error="std.error.T", statistic="statistic.T", p.value="p.value.T", lci="lci.T", uci="uci.T")
    
    f$estimate.T <- NULL
    f$std.error.T <- NULL
    f$statistic.T <- NULL
    f$p.value.T <- NULL
    f$lci.T <- NULL
    f$uci.T <- NULL
    f$subgroup <- F
    f <- f %>% dplyr::rename(estimate="estimate.F", std.error="std.error.F", statistic="statistic.F", p.value="p.value.F", lci="lci.F", uci="uci.F")

    d <- rbind(e, f)

    # create plot
    d$subgroup <- as.character(d$subgroup)
    d[which(d$subgroup == "FALSE")]$subgroup <- "Q1"
    d[which(d$subgroup == "TRUE")]$subgroup <- "Q2"
    p <- ggplot(d, aes(x=f, y=estimate, ymin=lci, ymax=uci, color=subgroup)) +
        coord_flip() +
        facet_grid(Trait~u, scales="free", space="free_y") +
        geom_point(size = 1.5) +
        geom_errorbar(width=.05) +
        theme_classic() +
        geom_rect(inherit.aes = F, show.legend = FALSE, data = tp, aes(fill = fill), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.15) +
        scale_y_continuous(breaks = scales::pretty_breaks(3)) +
        geom_hline(yintercept = c(0), linetype = "dashed", color = "grey") +
        scale_fill_manual(values=brewer.pal(2,"Paired")[1:2]) +
        scale_color_manual(values=brewer.pal(4,"Paired")[3:4], name="Subgroup") +
        theme(
            axis.title.y = element_blank(),
            strip.background = element_blank(),
            strip.text.y = element_text(angle = 0),
            legend.position = "bottom",
            panel.spacing.y = unit(0, "lines"),
            legend.box.background = element_rect(colour = "black")
        ) +
        ylab("Per allele effect estimate stratified by modifier, SD (95% CI)")
    return(p)
}

pdf("gxe-qual.pdf", width=8.5, height=8)
print(get_plot(get_dat("data/gxe-qual.txt")))
dev.off()