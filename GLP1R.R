library("ieugwasr")
library("data.table")
library("jlst")
source("funs.R")
set.seed(1234)
options(ieugwasr_api="http://64.227.44.193:8006/")

# load phenotypes
d <- fread("data/body_mass_index.21001.0.0.txt")

# extract cis-QTLs from opengwas
qtl <- ieugwasr::tophits(c("prot-a-1219", "finn-a-GLP1ANA", "eqtl-a-ENSG00000112164"))
qtl <- qtl %>% dplyr::filter(chr=="6")

# extract GTEx cis-eQTLs
gtex <- fread("data/GLP1R.csv")
gtex$ea <- str_split(gtex[['Variant Id']], "_", simplify = T)[,4]
gtex$oa <- str_split(gtex[['Variant Id']], "_", simplify = T)[,3]

# combine
df <- rbind(
    qtl %>% dplyr::select(beta, rsid, p, id, ea, nea) %>% dplyr::rename(pval="p", study="id", oa="nea") %>% dplyr::mutate(tissue="blood"),
    gtex %>% dplyr::select(NES, "SNP Id", "P-Value", ea, oa, Tissue) %>% dplyr::rename(pval="P-Value", beta="NES", rsid="SNP Id", tissue="Tissue") %>% dplyr::mutate(study="gtex")
)
df <- df %>% dplyr::filter(pval < 5e-8)

# clump
df <- ld_clump(df)

# test for vQTL effect
results <- data.frame()
for (i in 1:nrow(df)){
    variant <- bgen.load("/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr06.bgen",
        rsids=df$rsid[i]
    )
    dosage <- as.data.frame(
        apply(variant$data, 1, function(data) { return(data[,1]*0 + data[,2]*1 + data[,3]*2) })
    )
    dosage$appieu <- row.names(dosage)
    if (nrow(dosage) == 0){
        next
    }
    snpid <-paste0(names(dosage)[1], "_", variant$variants$allele0, "_", variant$variants$allele1)
    names(dosage)[1] <- snpid
    dosage <- merge(d, dosage, "appieu")
    p_main <- vartest(dosage %>% dplyr::pull("body_mass_index.21001.0.0"), dosage %>% dplyr::pull(!!snpid), covar = dosage %>% dplyr::select("age_at_recruitment.21022.0.0", "sex.31.0.0", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"), covar.var = T, type = 2, x.sq = T)$test$P
    p_log <- vartest(dosage %>% dplyr::pull("body_mass_index.21001.0.0") %>% log, dosage %>% dplyr::pull(!!snpid), covar = dosage %>% dplyr::select("age_at_recruitment.21022.0.0", "sex.31.0.0", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"), covar.var = T, type = 2, x.sq = T)$test$P
    results <- rbind(results, data.frame(
        p_main, p_log, rsid=snpid
    ))
}