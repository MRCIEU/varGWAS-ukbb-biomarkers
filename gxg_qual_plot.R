library("data.table")
library("ggplot2")
library("dplyr")
library("stringr")
library("ggpubr")
library("viridis")
library('forestplot')
library("viridis")
library("ggpubr")
library("RColorBrewer")
library("grid")
source("funs.R")
set.seed(123)

# read in gxg qual
d <- fread("data/gxg-qual.txt")
d$rsid.1 <- sapply(d$term, get_rsid)
d$rsid.2 <- sapply(d$mod_snp, get_rsid)
d$gene.1 <- sapply(d$rsid.1, get_gene)
d$gene.2 <- sapply(d$rsid.2, get_gene)
d$lci <- d$estimate - (1.96 * d$std.error)
d$uci <- d$estimate + (1.96 * d$std.error)
d$Trait <- sapply(d$trait, function(x) return(biomarkers_abr[biomarkers==x]))

d$f <- rep(NA, nrow(d))
for (i in 1:nrow(d)){
    if (is.na(d$gene.1[i])){
        d$Trait[i] <-  paste0(d$Trait[i], "\n", d$rsid.1[i], "\n",d$gene.2[i], " (", d$rsid.2[i], ")")
        d$f[i] <- paste0(d$rsid.1[i], " x ",d$gene.2[i], " (", d$rsid.2[i], ")")
    } else if (is.na(d$gene.2[i])){
        d$Trait[i] <-  paste0(d$Trait[i], "\n", d$gene.1[i], " (", d$rsid.1[i], ")\n", d$rsid.2[i])
        d$f[i] <- paste0(d$gene.1[i], " (", d$rsid.1[i], ") x ", d$rsid.2[i])
    } else{
        d$Trait[i] <-  paste0(d$Trait[i], "\n", d$gene.1[i], " (", d$rsid.1[i], ")\n",d$gene.2[i], " (", d$rsid.2[i], ")")
        d$f[i] <- paste0(d$gene.1[i], " (", d$rsid.1[i], ") x ",d$gene.2[i], " (", d$rsid.2[i], ")")
    }
}

d1 <- d %>% dplyr::select(mod, estimate, lci, uci, Trait, f) %>% dplyr::rename(copies="mod") %>% dplyr::mutate(int=F)
d1$model <- "Mean"

d <- fread("data/gxg-qual-var.txt")
d <- cbind(d, as.data.frame(str_split(d$term, ":", simplify=T), stringsAsFactors=F))
d$rsid.1 <- sapply(d$V1, get_rsid)
d$rsid.2 <- sapply(d$V2, get_rsid)
d$gene.1 <- sapply(d$rsid.1, get_gene)
d$gene.2 <- sapply(d$rsid.2, get_gene)
d$Trait <- sapply(d$trait, function(x) return(biomarkers_abr[biomarkers==x]))

d$f <- rep(NA, nrow(d))
for (i in 1:nrow(d)){
    if (is.na(d$gene.1[i])){
        d$Trait[i] <-  paste0(d$Trait[i], "\n", d$rsid.1[i], "\n",d$gene.2[i], " (", d$rsid.2[i], ")")
        d$f[i] <- paste0(d$rsid.1[i], " x ",d$gene.2[i], " (", d$rsid.2[i], ")")
    } else if (is.na(d$gene.2[i])){
        d$Trait[i] <-  paste0(d$Trait[i], "\n", d$gene.1[i], " (", d$rsid.1[i], ")\n", d$rsid.2[i])
        d$f[i] <- paste0(d$gene.1[i], " (", d$rsid.1[i], ") x ", d$rsid.2[i])
    } else{
        d$Trait[i] <-  paste0(d$Trait[i], "\n", d$gene.1[i], " (", d$rsid.1[i], ")\n",d$gene.2[i], " (", d$rsid.2[i], ")")
        d$f[i] <- paste0(d$gene.1[i], " (", d$rsid.1[i], ") x ",d$gene.2[i], " (", d$rsid.2[i], ")")
    }
}

d$phi_f <- NULL
d$phi_p <- NULL
e <- d
f <- d
e$estimate <- e$phi_x1
e$se <- e$se_x1
e$phi_x1 <- NULL
e$phi_x2 <- NULL
e$se_x1 <- NULL
e$se_x2 <- NULL
e$copies <- 1
f$estimate <- f$phi_x2
f$se <- f$se_x2
f$phi_x1 <- NULL
f$phi_x2 <- NULL
f$se_x1 <- NULL
f$se_x2 <- NULL
f$copies <- 2
d <- rbind(e,f)
d$lci <- d$estimate - (1.96 * d$se)
d$uci <- d$estimate + (1.96 * d$se)

d2 <- d %>% dplyr::select(copies, estimate, lci, uci, int, Trait, f)
d2$model <- "Variance"

# combine data
d <- rbind(d1, d2)
d$copies <- as.factor(d$copies)
d$int_t <- "Unadjusted"
d$int_t[d$int] <- "Adjusted"
d$int_t <- factor(d$int_t, levels=c("Unadjusted", "Adjusted"))

p <- ggplot(d, aes(x=copies, y=estimate, ymin=lci, ymax=uci, group=int_t, shape=int_t)) +
    geom_point(size = 2, position = position_dodge(width = 0.9)) +
    geom_errorbar(width=.05, position = position_dodge(width = 0.9)) +
    theme_classic() +
    facet_grid(Trait~model, scale="free") + 
    scale_y_continuous(breaks = scales::pretty_breaks(3)) +
    geom_hline(yintercept = c(0), linetype = "dashed", color = "grey") +
    theme(
        strip.background = element_blank(),
        panel.spacing.y = unit(1, "lines")
    ) +
    ylab("Trait (SD, 95% CI)") +
    xlab("SNP (genotype copies)") +
    labs(shape="Interaction")

pdf("gxg-qual.pdf", height=10)
print(p)
dev.off()