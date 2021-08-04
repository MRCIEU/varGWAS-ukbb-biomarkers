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

    # drop BMI
    d <- d %>% filter(trait != "body_mass_index.21001.0.0")

    # filter SNPs to show
    d <- d %>% filter(p.value < 5e-8)

    # merge
    d <- cbind(d, as.data.frame(str_split(d$term, ":", simplify=T), stringsAsFactors=F))
    d$y <- sapply(d$trait, function(x) return(biomarkers_abr[biomarkers==x]))
    d$u <- sapply(d$V1, get_trait_name)

    # map SNP to rsid & gene
    lookup <- fread("all.vqtls.txt")
    lookup$key <- paste0(lookup$key, "-", lookup$Trait)
    d$key <- paste0(d$V2, "-", d$y)
    d <- merge(d, lookup, "key")
    d$key <- NULL
    d$gene <- str_split(d[["Nearest Gene"]], ",", simplify=T)[,1]
    d$f <- paste0(d$gene, " (", d$RSID, d$EA, ")")
    d$u <- factor(d$u)
    levels(d$u) <- list(Age="Age At Recruitment", Sex="Sex", BMI="Body Mass Index", PA="Summed Minutes Activity", Alcohol="Alcohol Intake Frequency", Smoking="Smoking Status")

    # add key
    d$tt <- paste0(d$trait, ":", d$term)

    return(d)
}

get_plot <- function(d, leg_name, title_name){
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

    # sort data by locus name
    d$f <- factor(d$f, levels=unique(d$f) %>% sort(decreasing=T))

    # threshold rep P
    d$p_sens <- d$p_sens < 5e-5
    d$p_sens <- factor(d$p_sens, levels=c(T, F))

    # create plot
    p <- ggplot(d, aes(x=f, y=estimate, ymin=lci, ymax=uci, color=p_sens)) +
        coord_flip() +
        facet_grid(Trait~u, scales="free", space="free_y") +
        geom_point(size = 1.5) +
        geom_errorbar(width=.05) +
        geom_hline(yintercept = c(-0.05, 0, 0.05), linetype = "dashed", color = "grey") +
        theme_classic() +
        scale_y_continuous(limits = c(-.1, .1), breaks=c(-.1, 0, .1)) +
        geom_rect(inherit.aes = F, show.legend = FALSE, data = tp, aes(fill = fill), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.15) +
        scale_color_grey(name = paste0(leg_name, " P < 5 x 10-5")) +
        scale_fill_manual(values=brewer.pal(2,"Paired")) +
        ggtitle(title_name) +
        theme(
            axis.title.y = element_blank(),
            strip.background = element_blank(),
            strip.text.y = element_text(angle = 0),
            legend.position = "bottom",
            panel.spacing.y = unit(0, "lines")
        ) +
        ylab("Genotype (dosage) * modifier (SD) interaction effect estimate, SD (95% CI)")
    return(p)
}

# load gxe effects
additive <- get_dat("data/gxe.txt")
multiplicative <- get_dat("data/gxe-log.txt")

# append sensitivity P value
additive <- merge(additive, fread("data/gxe-log.txt") %>% mutate(tt=paste0(trait, ":", term)) %>% select(tt, p.value) %>% rename(p_sens="p.value"), "tt")
multiplicative <- merge(multiplicative, fread("data/gxe.txt") %>% mutate(tt=paste0(trait, ":", term)) %>% select(tt, p.value) %>% rename(p_sens="p.value"), "tt")

# save plot
pdf("gxe-additive.pdf", height=12, width=11)
print(get_plot(additive, "Multiplicative scale", "Gene-environment interaction effect (additive scale)"))
dev.off()

pdf("gxe-multiplicative.pdf", height=12, width=11)
print(get_plot(multiplicative, "Additive scale", "Gene-environment interaction effect (multiplicative scale)"))
dev.off()