library('data.table')
library('dplyr')
library('qqman')
library('stringr')
library('ggpubr')
source("funs.R")
set.seed(1234)

mean_man <- vector('list', length(biomarkers))
mean_qq <- vector('list', length(biomarkers))
var_man <- vector('list', length(biomarkers))
var_qq <- vector('list', length(biomarkers))
labs <- rep(NA, length(biomarkers))

for (i in 1:length(biomarkers)){
    trait_name <- biomarkers[i]
    trait_name <- str_split(trait_name, "\\.", simplify = TRUE)[,1]
    trait_name <- gsub("_", " ", trait_name)
    trait_name <- str_to_title(trait_name)
    labs[i] <- trait_name
    message(paste0("trait ", biomarkers[i]))
    message(paste0("trait name ", trait_name))

    # load vGWAS and QC
    data <- get_variants(biomarkers[i])

    # manhattan
    manhattan(data, ylim = c(0, 25), chr="chr", bp="pos", p="phi_p", snp="rsid")
    manhattan(data, ylim = c(0, 25), chr="chr", bp="pos", p="p", snp="rsid")

    # qq plot
    qq(data$phi_p)
    qq(data$p)
}

ggarrange(mean_man, labels = labs, ncol = 3, nrow = 10)