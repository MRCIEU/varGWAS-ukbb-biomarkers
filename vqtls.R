library("dplyr")
library("broom")
library("data.table")
library("stringr")
library("varGWASR")
source("funs.R")
set.seed(1234)
options(ieugwasr_api="http://64.227.44.193:8006/")

# load phenotypes
#load("data/pheno.RData")
#linker <- get_filtered_linker(drop_standard_excl=TRUE, drop_non_white_british=TRUE, drop_related=TRUE, application="15825")
#covariates <- get_covariates()
#pc <- get_genetic_principal_components()
#dat <- merge(linker, covariates, "appieu")
#dat <- merge(dat, pheno, by.x="app15825", by.y="eid")
#dat <- merge(dat, pc, "appieu")

# SD standardise traits
#for (trait in biomarkers){
#    dat[[trait]] <- dat[[trait]] / sd(dat[[trait]], na.rm=T)
#}

# load clumped vQTLs
#d <- fread("data/vqtls.txt")
#d$key <- paste0("chr", d$chr, "_", d$pos, "_", d$oa, "_", d$ea)
#d$chr_pos <- paste0(d$chr, ":", d$pos)

# load dosage
#snps <- d %>% dplyr::select(chr, pos, oa, ea)
#snps <- unique(snps)
#for (i in 1:nrow(snps)){
#    dosage <- tryCatch(
#        expr = {
#            extract_variant_from_bgen(as.character(snps$chr[i]), as.double(snps$pos[i]), snps$oa[i], snps$ea[i])
#        },
#        error = function(e){ 
#            NULL
#        }
#    )
#
#    if (is.null(dosage)){
#        warning(paste0("skipping variant: ", snps$chr[i], ":", snps$pos[i]))
#        next
#    }
#    dat <- merge(dat, dosage, "appieu")
#}

#save.image(file = "data/vqtls.RData")
load("data/vqtls.RData")

results <- data.frame()
for (i in 1:nrow(d)){
    message(paste0("working on: ", d$key[i]))
    tmp <- dat %>% dplyr::select(!!d$key[i], !!d$trait[i], c("age_at_recruitment.21022.0.0", "sex.31.0.0", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")) %>% tidyr::drop_na(.)
    tmp$PC1_xq <- tmp$PC1^2
    tmp$PC2_xq <- tmp$PC2^2
    tmp$PC3_xq <- tmp$PC3^2
    tmp$PC4_xq <- tmp$PC4^2
    tmp$PC5_xq <- tmp$PC5^2
    tmp$PC6_xq <- tmp$PC6^2
    tmp$PC7_xq <- tmp$PC7^2
    tmp$PC8_xq <- tmp$PC8^2
    tmp$PC9_xq <- tmp$PC9^2
    tmp$PC10_xq <- tmp$PC10^2
    tmp$age_at_recruitment.21022.0.0_xq <- tmp$age_at_recruitment.21022.0.0^2

    # drop Z > 5SD
    s <- sd(tmp[[d$trait[i]]])
    tmp$z <- (tmp[[d$trait[i]]] - mean (tmp[[d$trait[i]]])) / s
    tmp$z <- abs(tmp$z)
    tmp <- tmp %>% dplyr::filter(z < 5)
    tmp$z <- NULL
    
    # main
    test_main <- varGWASR::model(
        data=tmp,
        x=d$key[i], 
        y=d$trait[i], 
        covar1=c("age_at_recruitment.21022.0.0", "sex.31.0.0", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"),
        covar2=c("age_at_recruitment.21022.0.0", "sex.31.0.0", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")
    )
    test_main$snp=d$key[i]
    test_main$outcome=d$trait[i]

    results <- rbind(results, test_main)
}

# select fields for paper
results$key <- stringr::str_split(results$snp, "_", simplify=T) %>% 
    as.data.frame %>% 
    dplyr::mutate(V1=gsub("chr", "", V1)) %>% 
    dplyr::mutate(key=paste0(V1, ":",V2)) %>% 
    dplyr::pull(key)

# add RSID
rsid <- d %>% dplyr::select("key", "rsid")
results <- merge(results, rsid, by.x="snp", by.y="key")

# get data on nearest gene
ng <- fread("data/nearest.txt")
ng <- unique(ng)
ng$key <- paste0(ng$V1,":",ng$V2)
ng$V1 <- NULL
ng$V2 <- NULL
names(ng)[1] <- "gene"
ng <- ng %>%
  group_by_at(vars(key)) %>%
  summarize(gene = toString(gene)) %>%
  ungroup()
ng$gene <- gsub(", ", "|", ng$gene)

# results with nearest gene
all <- merge(results, ng, "key", all.x=T)
all$key <- NULL

# tidy outcome name
all$outcome <- sapply(all$outcome, function(x) biomarkers_abr[x==biomarkers], simplify=T)

# write to table
write.csv(all, file="Table S1.csv", quote=F, row.names=F)