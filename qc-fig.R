library('data.table')
library('dplyr')
library('qqman')
library('stringr')
source("funs.R")
set.seed(1234)

for (trait in c(
    "alanine_aminotransferase.30620.0.0",
    "albumin.30600.0.0",
    "alkaline_phosphatase.30610.0.0",
    "apolipoprotein_a.30630.0.0",
    "apolipoprotein_b.30640.0.0",
    "aspartate_aminotransferase.30650.0.0",
    "c_reactive_protein.30710.0.0",
    "calcium.30680.0.0",
    "cholesterol.30690.0.0",
    "creatinine.30700.0.0",
    "cystatin_c.30720.0.0",
    "direct_bilirubin.30660.0.0",
    "gamma_glutamyltransferase.30730.0.0",
    "glucose.30740.0.0",
    "glycated_haemoglobin.30750.0.0",
    "hdl_cholesterol.30760.0.0",
    "igf_1.30770.0.0",
    "ldl_direct.30780.0.0",
    "lipoprotein_a.30790.0.0",
    "oestradiol.30800.0.0",
    "phosphate.30810.0.0",
    "rheumatoid_factor.30820.0.0",
    "shbg.30830.0.0",
    "testosterone.30850.0.0",
    "total_bilirubin.30840.0.0",
    "total_protein.30860.0.0",
    "triglycerides.30870.0.0",
    "urate.30880.0.0",
    "urea.30670.0.0",
    "vitamin_d.30890.0.0"
)){
    trait_name <- trait
    trait_name <- str_split(trait_name, "\\.", simplify = TRUE)[,1]
    trait_name <- gsub("_", " ", trait_name)
    trait_name <- str_to_title(trait_name)
    message(paste0("trait ", trait))
    message(paste0("trait name ", trait_name))

    # load vGWAS and QC
    data <- get_variants(trait)

    # manhattan
    png(paste0("data/", trait, "_phi_manhattan.png"))
    manhattan(data, ylim = c(0, 25), chr="chr", bp="pos", p="phi_p", snp="rsid", main = trait_name)
    dev.off()
    png(paste0("data/", trait, "_beta_manhattan.png"))
    manhattan(data, ylim = c(0, 25), chr="chr", bp="pos", p="p", snp="rsid", main = trait_name)
    dev.off()

    # qq plot
    png(paste0("data/", trait, "_phi_qq.png"))
    qq(data$phi_p, main = trait_name)
    dev.off()
    png(paste0("data/", trait, "_beta_qq.png"))
    qq(data$p, main = trait_name)
    dev.off()
}