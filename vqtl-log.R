library("data.table")
library("varGWASR")
source("funs.R")
set.seed(1234)

covar <- c("sex.31.0.0", "age_at_recruitment.21022.0.0", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")

# load clumped vQTLs
d <- fread("data/vqtls.txt")

# replicate on log scale
results <- data.frame()
for (trait in unique(d$trait)){

    # subset trait
    e <- d %>% dplyr::filter(trait == !!trait)

    # load data for trait
    f <- fread(paste0("data/", trait, ".txt"))
    
    # add in log scale
    f[[paste0(trait, "_log")]] <- log(f[[trait]])

    # loop over vQTLs and test for effect
    for (i in 1:nrow(e)){
        # load vQTL
        dosage <- extract_variant_from_bgen(as.character(e$chr[i]), as.double(e$pos[i]), e$oa[i], e$ea[i])
        g <- merge(f, dosage, "appieu")
        g <- na.omit(g)

        # test effect
        snp <- paste0("chr", e$chr[i], "_", e$pos[i], "_", e$oa[i], "_", e$ea[i])
        fit_main <- varGWASR::model(g, snp, trait, covar1=covar, covar2=covar)
        names(fit_main) <- paste0(names(fit_main), ".main")
        fit_log <- varGWASR::model(g, snp, paste0(trait, "_log"), covar1=covar, covar2=covar)
        names(fit_log) <- paste0(names(fit_log), ".log")
        fit <- cbind(fit_main, fit_log)
        fit$trait <- trait
        fit$snp <- snp

        # save
        results <- rbind(results, fit)
    }

}

write.csv(results, file="data/vqtl-log.csv")