library("data.table")
library("stringr")
source("funs.R")
set.seed(13)

map_snps <- function(row){
    snps <- grep("chr", stringr::str_split(row, " ")[[1]], value=T)
    return(paste0(snps, collapse="|"))
}

# gxg
gxg <- fread("data/gxg-add-finemap.txt", select=c("term", "trait", "formula"))
gxg$snps <- sapply(gxg$formula, function(row) map_snps(row))
gxg$formula <- NULL

# gxe
gxe <- fread("data/gxe-finemap.txt", select=c("term", "trait", "formula"))
gxe$snps <- sapply(gxe$formula, function(row) map_snps(row))
gxe$formula <- NULL

# all
all <- rbind(gxg, gxe)
all$snp <- stringr::str_split(all$term, ":", simplify=T)[,2]
all$modifier <- stringr::str_split(all$term, ":", simplify=T)[,1]
all$modifier <- gsub("age_at_recruitment.21022.0.0", "Age", all$modifier)
all$modifier <- gsub("alcohol_intake_frequency.1558.0.0", "Alcohol", all$modifier)
all$modifier <- gsub("body_mass_index.21001.0.0", "BMI", all$modifier)
all$modifier <- gsub("sex.31.0.0", "Sex", all$modifier)
all$modifier <- gsub("smoking_status.20116.0.0", "Smoking", all$modifier)
all$term <- NULL
lookup <- fread("Table S1.csv", select=c("snp", "gene", "rsid"))
lookup <- unique(lookup)
all <- merge(all, lookup,"snp")
all$gene <- stringr::str_split(all$gene, "\\|", simplify=T)[,1]
all$y <- sapply(all$trait, function(z) biomarkers_abr[z == biomarkers])

all$trait <- NULL
all$x <- paste0(all$gene, "(", all$rsid, ")")
all$snp <- NULL

all <- all %>% dplyr::select(x, modifier, y, snps)

names(all) <- c(
    "snp",
    "modifier",
    "outcome",
    "fine_mapped"
)

write.csv(all, file="Table S3.csv", quote=F, row.names=F)