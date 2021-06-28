library("data.table")
library("ggplot2")
library("dplyr")
library("stringr")
library("ggpubr")
library("viridis")
library('forestplot')
library("RColorBrewer")
library("grid")
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
wf$key <- paste0(wf$term, "-", wf$trait)
wf$trait <- NULL
wf$term <- NULL

# filter genome-wide sig hits
d <- d %>% filter(trait != "body_mass_index.21001.0.0")
d <- d %>% filter(p.value < 5e-8)
d$key <- paste0(d$term, "-", d$trait)

# merge
d <- merge(d, wf, "key")
d$key <- NULL
d <- cbind(d, as.data.frame(str_split(d$term, ":", simplify=T), stringsAsFactors=F))
d$y <- sapply(d$trait, function(x) return(biomarkers_abr[biomarkers==x]))
d$u <- sapply(d$V1, get_trait_name)

# map SNP to rsid & gene
lookup <- fread("all.vqtls.txt")
lookup$key <- paste0(lookup$key, "-", lookup$Trait)
d$key <- paste0(d$V2, "-", d$y)
d <- merge(d, lookup, "key")
d$key <- NULL
d$f <- as.factor(d$RSID)
d$u <- factor(d$u)
levels(d$u) <- list(Age="Age At Recruitment", Smoking="Smoking Status", Sex="Sex", BMI="Body Mass Index", Alcohol="Alcohol Intake Frequency")

# select fields
main <- d %>% select(Trait, f, estimate, lci.x, uci.x, u) %>% rename(lci="lci.x", uci="uci.x")
main$analysis <- "Main"
wf <- d %>% select(Trait, f, u, beta, lci.y, uci.y) %>% rename(estimate="beta", lci="lci.y", uci="uci.y")
wf$analysis <- "WF"
d <- rbind(main, wf)

# plot
get_plot <- function(data){
    p <- ggplot(data=data, aes(x=f, y=estimate, ymin=lci, ymax=uci, group=u, color=u)) +
        coord_flip() +
        facet_wrap(Trait~., scales = "free")+
        geom_point(position=position_dodge(width=1)) +
        geom_errorbar(width=.05, position=position_dodge(width=1)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
        theme_classic() +
        xlab("SNP") + 
        ylab("Genotype (dosage) * modifier (SD) interaction effect estimate, SD (95% CI)") +
        scale_y_continuous(breaks = scales::breaks_pretty(n=4)) +
        theme(
            strip.text.y = element_text(angle=0),
            strip.background = element_blank(),
            axis.title.y = element_blank()
        ) +
        labs(color = "Modifier")
        return(p)
}

pdf("gxe-main.pdf", width=12, height=10)
print(get_plot(main))
dev.off()

pdf("gxe-wf.pdf", width=12, height=10)
print(get_plot(wf))
dev.off()