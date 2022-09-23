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
d <- fread("data/gxe-qual-var1.txt")
d$mod_pheno <- gsub("body_mass_index.21001.0.0_b", "BMI", d$mod_pheno)
d$mod_pheno <- gsub("sex.31.0.0_b", "Sex", d$mod_pheno)
d$rsid.1 <- sapply(d$term, get_rsid)
d$gene.1 <- sapply(d$rsid.1, get_gene)
d$lci <- d$estimate - (1.96 * d$std.error)
d$uci <- d$estimate + (1.96 * d$std.error)
d$Trait <- sapply(d$trait, function(x) return(biomarkers_abr[biomarkers==x]))
d$Trait <-  paste0(d$Trait, "\n", d$gene.1, " (", d$rsid.1, ")\n",d$mod_pheno)

d1 <- d %>% dplyr::select(mod, estimate, lci, uci, Trait) %>% dplyr::rename(copies="mod") %>% dplyr::mutate(int=F)
d1$model <- "Mean"
d1$copies[grepl("BMI", d1$Trait) & d1$copies ==0] <- "Low BMI"
d1$copies[grepl("BMI", d1$Trait) & d1$copies ==1] <- "High BMI"
d1$copies[grepl("Sex", d1$Trait) & d1$copies ==0] <- "Female"
d1$copies[grepl("Sex", d1$Trait) & d1$copies ==1] <- "Male"

d <- fread("data/gxe-qual-var2.txt")
d$phi_f <- NULL
d$phi_p <- NULL
e <- d
f <- d
e$phi <- e$phi_x1
e$se <- e$se_x1
e$phi_x1 <- NULL
e$phi_x2 <- NULL
e$se_x1 <- NULL
e$se_x2 <- NULL
e$gt <- 1
f$phi <- f$phi_x2
f$se <- f$se_x2
f$phi_x1 <- NULL
f$phi_x2 <- NULL
f$se_x1 <- NULL
f$se_x2 <- NULL
f$gt <- 2
d <- rbind(e,f)

d <- cbind(d, as.data.frame(str_split(d$term, ":", simplify=T), stringsAsFactors=F))
d$V2 <- gsub("body_mass_index.21001.0.0", "BMI", d$V2)
d$V2 <- gsub("sex.31.0.0", "Sex", d$V2)
d$rsid.1 <- sapply(d$V1, get_rsid)
d$gene.1 <- sapply(d$rsid.1, get_gene)
d$lci <- d$phi - (1.96 * d$se)
d$uci <- d$phi + (1.96 * d$se)
d$Trait <- sapply(d$trait, function(x) return(biomarkers_abr[biomarkers==x]))
d$Trait <-  paste0(d$Trait, "\n", d$gene.1, " (", d$rsid.1, ")\n",d$V2)
d$gt <- as.factor(d$gt)

d2 <- d %>% dplyr::select(gt, phi, lci, uci, int, Trait) %>% dplyr::rename(copies="gt",estimate="phi")
d2$model <- "Variance"

# combine data
d <- rbind(d1, d2)
d$int_t <- "Unadjusted"
d$int_t[d$int] <- "Adjusted"
d$int_t <- factor(d$int_t, levels=c("Unadjusted", "Adjusted"))

p <- ggplot(d, aes(x=copies, y=estimate, ymin=lci, ymax=uci, group=int_t, shape=int_t)) +
    geom_point(size = 2, position = position_dodge(width = 0.9)) +
    geom_errorbar(width=.05, position = position_dodge(width = 0.9)) +
    theme_classic() +
    facet_wrap(Trait~model, ncol=2, scales="free") + 
    scale_y_continuous(breaks = scales::pretty_breaks(3)) +
    geom_hline(yintercept = c(0), linetype = "dashed", color = "grey") +
    theme(
        strip.background = element_blank(),
        legend.background = element_blank(),
        axis.title.x=element_blank()
    ) +
    ylab("Trait (SD, 95% CI)") +
    labs(shape="Interaction")

pdf("gxe-qual-var.pdf", height=10)
print(p)
dev.off()