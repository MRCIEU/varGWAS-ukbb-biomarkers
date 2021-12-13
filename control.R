load("data/pheno.RData")

library('data.table')
library('dplyr')
library("varGWASR")
library("ggplot2")
source("funs.R")
set.seed(1234)

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

# SD scale outcomes
dat$standing_height.50.0.0 <- dat$standing_height.50.0.0 / sd(dat$standing_height.50.0.0, na.rm=T)
dat$forced_expiratory_volume_best_measure.20150.0.0 <- dat$forced_expiratory_volume_best_measure.20150.0.0 / sd(dat$forced_expiratory_volume_best_measure.20150.0.0, na.rm=T)
dat$forced_vital_capacity_best_measure.20151.0.0 <- dat$forced_vital_capacity_best_measure.20151.0.0 / sd(dat$forced_vital_capacity_best_measure.20151.0.0, na.rm=T)

# extract smoking heaviness SNP
dosage <- extract_variant_from_bgen("15", 78894339, "G", "A")
dat <- merge(dat, dosage, "appieu")
dat$int <- dat$chr15_78894339_G_A * dat$smoking_status.20116.0.0

# test for effect of SNP on trait variance w/wo adjusting for SNP * smoking status
covar <- c("age_at_recruitment.21022.0.0", "sex.31.0.0", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")
results <- data.frame()

fitfun <- function(d, x, y, covar, trait, interaction){
    fit <- varGWASR::model(d, x, y, covar1 = covar, covar2 = covar)
    res <- data.frame()
    res <- rbind(res, data.frame(
        phi=fit$phi_x1,
        lci=fit$phi_x1 - (1.96 * fit$se_x1),
        uci=fit$phi_x1 + (1.96 * fit$se_x1),
        gt=1,
        trait,
        interaction
    ))
    res <- rbind(res, data.frame(
        phi=fit$phi_x2,
        lci=fit$phi_x2 - (1.96 * fit$se_x2),
        uci=fit$phi_x2 + (1.96 * fit$se_x2),
        gt=2,
        trait,
        interaction
    ))
    return(res)
}

## FVC
fvc <- dat %>% dplyr::select(chr15_78894339_G_A, int, forced_vital_capacity_best_measure.20151.0.0, smoking_status.20116.0.0, age_at_recruitment.21022.0.0, sex.31.0.0, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10) %>% na.omit
results <- rbind(results, fitfun(fvc, "chr15_78894339_G_A", "forced_vital_capacity_best_measure.20151.0.0", covar, "FVC", "Unadjusted"))
results <- rbind(results, fitfun(fvc, "chr15_78894339_G_A", "forced_vital_capacity_best_measure.20151.0.0", c(covar, "int", "smoking_status.20116.0.0"), "FVC", "Adjusted"))

## FEV1
fev1 <- dat %>% dplyr::select(chr15_78894339_G_A, int, forced_expiratory_volume_best_measure.20150.0.0, smoking_status.20116.0.0, age_at_recruitment.21022.0.0, sex.31.0.0, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10) %>% na.omit
results <- rbind(results, fitfun(fev1, "chr15_78894339_G_A", "forced_expiratory_volume_best_measure.20150.0.0", covar, "FEV1", "Unadjusted"))
results <- rbind(results, fitfun(fev1, "chr15_78894339_G_A", "forced_expiratory_volume_best_measure.20150.0.0", c(covar, "int", "smoking_status.20116.0.0"), "FEV1", "Adjusted"))

## Height
height <- dat %>% dplyr::select(chr15_78894339_G_A, int, standing_height.50.0.0, smoking_status.20116.0.0, age_at_recruitment.21022.0.0, sex.31.0.0, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10) %>% na.omit
results <- rbind(results, fitfun(height, "chr15_78894339_G_A", "standing_height.50.0.0", covar, "Height", "Unadjusted"))
results <- rbind(results, fitfun(height, "chr15_78894339_G_A", "standing_height.50.0.0", c(covar, "int", "smoking_status.20116.0.0"), "Height", "Adjusted"))

# plot
results$gt <- as.factor(results$gt)
results$trait <- factor(results$trait, levels = c("FEV1", "FVC", "Height"))

pdf("control.pdf")
ggplot(results, aes(x=gt, y=phi, ymin=lci, ymax=uci, group=gt)) +
    geom_point(size = 1.5, position = position_dodge(width = 0.9)) +
    geom_errorbar(width=.05, position = position_dodge(width = 0.9)) +
    theme_classic() +
    facet_grid(trait~interaction, scale="free") + 
    scale_y_continuous(breaks = scales::pretty_breaks(3)) +
    geom_hline(yintercept = c(0), linetype = "dashed", color = "grey") +
    theme(
        strip.background = element_blank(),
        panel.spacing.y = unit(1, "lines")
    ) +
    ylab("Trait variance (95% CI)") +
    xlab("SNP (genotype copies)")
dev.off()