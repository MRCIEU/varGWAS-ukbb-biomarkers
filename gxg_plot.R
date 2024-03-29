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
options(ieugwasr_api="http://web-dc1-bms-d0.infra.bris.ac.uk:5002/")

get_dat <- function(file){
    # read in gxe results
    d <- fread(file)
    d$lci <- d$estimate - (d$std.error * 1.96)
    d$uci <- d$estimate + (d$std.error * 1.96)

    # drop BMI
    d <- d %>% dplyr::filter(trait != "body_mass_index.21001.0.0")

    # filter SNPs to show
    d <- d %>% dplyr::filter(p.value < 5e-8)

    # merge
    d <- cbind(d, as.data.frame(str_split(d$term, ":", simplify=T), stringsAsFactors=F))
    d$Trait <- sapply(d$trait, function(x) return(biomarkers_abr[biomarkers==x]))
    d$rsid.1 <- sapply(d$V1, get_rsid)
    d$rsid.2 <- sapply(d$V2, get_rsid)
    d$gene.1 <- sapply(d$rsid.1, get_gene)
    d$gene.2 <- sapply(d$rsid.2, get_gene)

    # add key
    d$tt <- paste0(d$trait, ":", d$term)
    d$f <- rep(NA, nrow(d))
    for (i in 1:nrow(d)){
        if (is.na(d$gene.1[i])){
            d$f[i] <- paste0(d$rsid.1[i], " x ",d$gene.2[i], " (", d$rsid.2[i], ")")
        } else if (is.na(d$gene.2[i])){
            d$f[i] <- paste0(d$gene.1[i], " (", d$rsid.1[i], ") x ", d$rsid.2[i])
        } else{
            d$f[i] <- paste0(d$gene.1[i], " (", d$rsid.1[i], ") x ",d$gene.2[i], " (", d$rsid.2[i], ")")
        }
    }

    return(d)
}

get_plot <- function(d, leg_name){
    # sort data by locus name
    d$f <- factor(d$f, levels=unique(d$f) %>% sort(decreasing=T))
    
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
            Trait=tr,
            fill=which(tr == unique(d$Trait)) %% 2
        ))
    }
    tp$fill <- as.factor(tp$fill)

    # sort data by locus name
    d$f <- factor(d$f, levels=unique(d$f) %>% sort(decreasing=F))

    # threshold rep P
    d$p_sens <- d$p_sens < 5e-8
    d$p_sens <- factor(d$p_sens, levels=c(T, F))

    # create plot
    p <- ggplot(d, aes(x=f, y=estimate, ymin=lci, ymax=uci, color=p_sens)) +
        coord_flip() +
        facet_grid(Trait~., scales="free", space="free_y") +
        geom_point(size = 1.5) +
        geom_errorbar(width=.05) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
        theme_classic() +
        scale_y_continuous(limits = c(-.1, .1), breaks=c(-.1, -0.05, 0, 0.05, .1)) +
        geom_rect(inherit.aes = F, show.legend = FALSE, data = tp, aes(fill = fill), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.15) +
        scale_color_grey(name = leg_name) +
        scale_fill_manual(values=brewer.pal(2,"Paired")) +
        theme(
            axis.title.y = element_blank(),
            strip.background = element_blank(),
            strip.text.y = element_text(angle = 0),
            legend.position = "bottom",
            legend.box.background = element_rect(colour = "black"),
            panel.spacing.y = unit(0, "lines")
        ) +
        ylab("Genotype * genotype (dosage) interaction effect estimate, SD (95% CI)")
    return(p)
}

# load gxe effects
additive <- get_dat("data/gxg.txt")
finemapped <- get_dat("data/gxg-add-finemap.txt")
multiplicative <- get_dat("data/gxg-log.txt")

# append sensitivity P value
additive <- merge(additive, fread("data/gxg-log.txt") %>% mutate(tt=paste0(trait, ":", term)) %>% select(tt, p.value) %>% rename(p_sens="p.value"), "tt")
multiplicative <- merge(multiplicative, fread("data/gxg.txt") %>% mutate(tt=paste0(trait, ":", term)) %>% select(tt, p.value) %>% rename(p_sens="p.value"), "tt")
finemapped$p_sens=1

# save plot
pdf("gxg-additive.pdf", height=6, width=8)
print(get_plot(additive, "Multiplicative (P < 5e-8)"))
dev.off()

pdf("gxg-multiplicative.pdf", height=5, width=8)
print(get_plot(multiplicative, "Additive (P < 5e-8)"))
dev.off()

pdf("gxg-finemap.pdf", height=6, width=8)
print(get_plot(finemapped, ""))
dev.off()