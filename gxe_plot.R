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

get_dat <- function(file, threshold=5e-8){
    # read in gxe results
    d <- fread(file)
    d$lci <- d$estimate - (d$std.error * 1.96)
    d$uci <- d$estimate + (d$std.error * 1.96)

    # drop BMI
    d <- d %>% dplyr::filter(trait != "body_mass_index.21001.0.0")

    # filter SNPs to show
    d <- d %>% dplyr::filter(p.value < threshold)

    # merge
    d <- cbind(d, as.data.frame(str_split(d$term, ":", simplify=T), stringsAsFactors=F))
    d$y <- sapply(d$trait, function(x) return(biomarkers_abr[biomarkers==x]))
    d$u <- sapply(d$V1, get_trait_name)

    # map SNP to rsid & gene
    lookup <- fread("Table S1.csv", select=c("snp", "gene", "rsid"))
    lookup <- unique(lookup)
    d <- merge(d, lookup, by.x="V2", by.y="snp")
    d$gene <- stringr::str_split(d$gene, "\\|", simplify=T)[,1]
    d$u <- factor(d$u)
    d$f <- paste0(d$gene, " (", d$rsid, ")")
    
    levels(d$u) <- list(Age="Age At Recruitment", Sex="Sex", BMI="Body Mass Index", PA="Summed Minutes Activity", Alcohol="Alcohol Intake Frequency", Smoking="Smoking Status")

    # add key
    d$tt <- paste0(d$trait, ":", d$term)
    d$Trait <- d$y

    return(d)
}

get_plot <- function(d, leg_name){
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
    d$p_sens <- d$p_sens < 5e-8
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
        scale_color_grey(name = leg_name) +
        scale_fill_manual(values=brewer.pal(2,"Paired")) +
        theme(
            axis.title.y = element_blank(),
            strip.background = element_blank(),
            strip.text.y = element_text(angle = 0),
            legend.position = "bottom",
            legend.background = element_blank(),
            legend.box.background = element_rect(colour = "black"),
            panel.spacing.y = unit(0, "lines")
        ) +
        ylab("Genotype (dosage) * modifier (SD) interaction effect estimate, SD (95% CI)")
    return(p)
}

# load gxe effects
additive <- get_dat("data/gxe.txt")
multiplicative <- get_dat("data/gxe-log.txt")
finemapped <- get_dat("data/gxe-finemap.txt", threshold=1)
finemapped$p_sens <- FALSE

# append sensitivity P value
additive <- merge(additive, fread("data/gxe-log.txt") %>% mutate(tt=paste0(trait, ":", term)) %>% select(tt, p.value) %>% rename(p_sens="p.value"), "tt")
multiplicative <- merge(multiplicative, fread("data/gxe.txt") %>% mutate(tt=paste0(trait, ":", term)) %>% select(tt, p.value) %>% rename(p_sens="p.value"), "tt")

# save plot
pdf("gxe-additive.pdf", height=12, width=11)
print(get_plot(additive, "Multiplicative (P < 5e-8)"))
dev.off()

pdf("gxe-multiplicative.pdf", height=12, width=11)
print(get_plot(multiplicative, "Additive (P < 5e-8)"))
dev.off()

pdf("gxe-finemapped.pdf", height=12, width=11)
print(get_plot(finemapped, ""))
dev.off()

# GxE table
additive$ea <- stringr::str_split(additive$V2, "_", simplify=T)[,4]
additive_tbl <- additive %>% 
    filter(tt != "hdl_cholesterol.30760.0.0:sex.31.0.0:chr19_45413576_C_T") %>%
    filter(tt != "triglycerides.30870.0.0:body_mass_index.21001.0.0:chr22_44324727_C_G") %>%
    filter(tt != "urate.30880.0.0:sex.31.0.0:chr4_10402838_T_C") %>%
    filter(p_sens < 5e-8) %>%
    arrange(p.value) %>%
    select(rsid, ea, u, y, estimate, lci, uci, p.value, gene)
additive_tbl$SNP <- paste0(additive_tbl$rsid, additive_tbl$ea)
additive_tbl$rsid <- NULL
additive_tbl$ea <- NULL
names(additive_tbl) <- c("Modifier", "Outcome", "Effect", "L95CI", "U95CI", "P", "Gene", "SNP")
additive_tbl <- additive_tbl %>% select(c("SNP", "Modifier", "Outcome", "Effect", "L95CI", "U95CI", "P", "Gene"))
additive_tbl$Effect <- round(additive_tbl$Effect, 2)
additive_tbl$L95CI <- round(additive_tbl$L95CI, 2)
additive_tbl$U95CI <- round(additive_tbl$U95CI, 2)
additive_tbl$P <- signif(additive_tbl$P, digits=2)
write.csv(additive_tbl, file="Table_2.csv", row.names=F, quote=F)