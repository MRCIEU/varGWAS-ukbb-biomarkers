load("data/pheno2.RData")
library('optparse')
library('data.table')
library('dplyr')
library('broom')
library('ieugwasr')
library("ggpubr")
library("lmtest")
library("ggplot2")
library("sandwich")
source("funs.R")
set.seed(1234)
options(ieugwasr_api="http://64.227.44.193:8006/")

# load linker
linker <- get_filtered_linker(drop_standard_excl=TRUE, drop_non_white_british=TRUE, drop_related=TRUE, application="15825")

# load covariates
covariates <- get_covariates()
covariates$chip <- as.numeric(as.factor(covariates$chip)) - 1
pc <- get_genetic_principal_components()

# merge data
dat <- merge(linker, covariates, "appieu")
dat <- merge(dat, pheno, by.x="app15825", by.y="eid")
dat <- merge(dat, pc, "appieu")

# split env into binary groups
dat$body_mass_index.21001.0.0_b <- dat$body_mass_index.21001.0.0 > median(dat$body_mass_index.21001.0.0, na.rm=T)
dat$sex.31.0.0_b <- dat$sex.31.0.0 > 0

# SD outcome
dat$hdl_cholesterol.30760.0.0 <- dat$hdl_cholesterol.30760.0.0 / sd(dat$hdl_cholesterol.30760.0.0, na.rm=T)
dat$triglycerides.30870.0.0 <- dat$triglycerides.30870.0.0 / sd(dat$triglycerides.30870.0.0, na.rm=T)
dat$urate.30880.0.0 <- dat$urate.30880.0.0 / sd(dat$urate.30880.0.0, na.rm=T)

# load dosages
dosage <- extract_variant_from_bgen("19", 45413576, "C", "T")
dat <- merge(dat, dosage, "appieu")
dosage <- extract_variant_from_bgen("22", 44324727, "C", "G")
dat <- merge(dat, dosage, "appieu")
dosage <- extract_variant_from_bgen("4", 10402838, "T", "C")
dat <- merge(dat, dosage, "appieu")

results <- data.frame()

get_lm_estimate <- function(f, k, dat, trait, g1, g2){
    mod <- lm(f, data=dat %>% dplyr::filter(!!sym(k) == T))
    t <- coeftest(mod, vcov = vcovHC(mod, type = "HC0")) %>% tidy
    res1 <- t[2,]
    mod <- lm(f, data=dat %>% dplyr::filter(!!sym(k) == F))
    t <- coeftest(mod, vcov = vcovHC(mod, type = "HC0")) %>% tidy
    res2 <- t[2,]
    res1$group <- g1
    res2$group <- g2
    res <- rbind(res1, res2)
    res$trait <- trait
    res$binary = FALSE
    res$k <- k
    return(res)
}

get_glm_estimate <- function(f, k, dat, trait, g1, g2){
    mod <- glm(f, family="binomial", data=dat %>% dplyr::filter(!!sym(k) == T))
    t <- tidy(mod)
    res1 <- t[2,]
    mod <- glm(f, family="binomial", data=dat %>% dplyr::filter(!!sym(k) == F))
    t <- tidy(mod)
    res2 <- t[2,]
    res1$group <- g1
    res2$group <- g2
    res <- rbind(res1, res2)
    res$trait <- trait
    res$binary = TRUE
    res$k <- k
    return(res)
}

# Effect of rs75627662 on HDL/CVD by sex
f <- as.formula(hdl_cholesterol.30760.0.0 ~ chr19_45413576_C_T + age_at_recruitment.21022.0.0 + sex.31.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10)
results <- rbind(results, get_lm_estimate(f, "sex.31.0.0_b", dat, "hdl_cholesterol.30760.0.0", "Male", "Female"))
f <- as.formula(vascular_problems.6150 ~ chr19_45413576_C_T + age_at_recruitment.21022.0.0 + sex.31.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10)
results <- rbind(results, get_glm_estimate(f, "sex.31.0.0_b", dat, "vascular_problems.6150", "Male", "Female"))

# Effect of rs738409 on TG/CVD by BMI
f <- as.formula(triglycerides.30870.0.0 ~ chr22_44324727_C_G + age_at_recruitment.21022.0.0 + sex.31.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10)
results <- rbind(results, get_lm_estimate(f, "body_mass_index.21001.0.0_b", dat, "triglycerides.30870.0.0", "BMI > 26.7 kg/m2", "BMI < 26.7 kg/m2"))

f <- as.formula(vascular_problems.6150 ~ chr22_44324727_C_G + age_at_recruitment.21022.0.0 + sex.31.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10)
results <- rbind(results, get_glm_estimate(f, "body_mass_index.21001.0.0_b", dat, "vascular_problems.6150", "BMI > 26.7 kg/m2", "BMI < 26.7 kg/m2"))

# Effect of rs4530622 on Urate/Gout by sex
f <- as.formula(urate.30880.0.0 ~ chr4_10402838_T_C + age_at_recruitment.21022.0.0 + sex.31.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10)
results <- rbind(results, get_lm_estimate(f, "sex.31.0.0_b", dat, "urate.30880.0.0", "Male", "Female"))

f <- as.formula(gout ~ chr4_10402838_T_C + age_at_recruitment.21022.0.0 + sex.31.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10)
results <- rbind(results, get_glm_estimate(f, "sex.31.0.0_b", dat, "gout", "Male", "Female"))

# plots
results$lci <- results$estimate - (1.96 * results$std.error)
results$uci <- results$estimate + (1.96 * results$std.error)
results$estimate[results$binary] <- exp(results$estimate[results$binary])
results$lci[results$binary] <- exp(results$lci[results$binary])
results$uci[results$binary] <- exp(results$uci[results$binary])

get_lm_plot <- function(results, x, y, title){
    results$group <- as.factor(results$group)
    p <- ggplot(results, aes(x=group, y=estimate, ymin=lci, ymax=uci)) +
        geom_point(size = 1.5, position = position_dodge(width = 0.9)) +
        geom_errorbar(width=.05, position = position_dodge(width = 0.9)) +
        theme_classic() +
        scale_y_continuous(breaks = scales::pretty_breaks(5)) +
        geom_hline(yintercept = c(0), linetype = "dashed", color = "grey") +
        theme(
            panel.spacing.y = unit(0, "lines"),
            legend.box.background = element_rect(colour = "black"),
            legend.position = "none"
        ) +
        ylab(x) +
        xlab(y) + 
        ggtitle(title)
        return(p)
}

get_glm_plot <- function(results, x, y, title){
    results$group <- as.factor(results$group)
    p <- ggplot(results, aes(x=group, y=estimate, ymin=lci, ymax=uci)) +
        geom_point(size = 1.5, position = position_dodge(width = 0.9)) +
        geom_errorbar(width=.05, position = position_dodge(width = 0.9)) +
        theme_classic() +
        scale_y_continuous(breaks = scales::pretty_breaks(5)) +
        geom_hline(yintercept = c(1), linetype = "dashed", color = "grey") +
        theme(
            panel.spacing.y = unit(0, "lines"),
            legend.box.background = element_rect(colour = "black"),
            legend.position = "none"
        ) +
        ylab(x) +
        xlab(y) + 
        ggtitle(title)
        return(p)
}

p1 <- get_lm_plot(results %>% filter(term == "chr19_45413576_C_T" & !binary), "HDL (95% CI)", "Sex", "APOE (rs75627662)")
p2 <- get_glm_plot(results %>% filter(term == "chr19_45413576_C_T" & binary), "Vascular disease (OR, 95% CI)", "Sex", "")

p3 <- get_lm_plot(results %>% filter(term == "chr22_44324727_C_G" & !binary), "TG (95% CI)", "BMI", "PNPLA3 (rs738409)")
p4 <- get_glm_plot(results %>% filter(term == "chr22_44324727_C_G" & binary), "Vascular disease (OR, 95% CI)", "BMI", "")

p5 <- get_lm_plot(results %>% filter(term == "chr4_10402838_T_C" & !binary), "Urate (95% CI)", "Sex", "SLC2A9 (rs4530622)")
p6 <- get_glm_plot(results %>% filter(term == "chr4_10402838_T_C" & binary), "Gout (OR, 95% CI)", "Sex", "")

p <- ggarrange(p1, p2, p3, p4, p5, p6, labels = c("A", "B", "C", "D", "E", "F"), ncol = 2, nrow = 3)

pdf("gxe-disease-qual.pdf", height=7*(2/3)*3, width=14)
print(p)
dev.off()