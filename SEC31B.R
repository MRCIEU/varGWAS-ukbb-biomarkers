library('data.table')
library('dplyr')
library('ieugwasr')
library('broom')
library('TwoSampleMR')
source("funs.R")
set.seed(124)

### MR effect of AST on SEC31B expression ###
# Get instruments
exposure_dat <- extract_instruments("ukb-d-30650_irnt")

# Get effects of instruments on outcome
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="eqtl-a-ENSG00000075826")

# Harmonise the exposure and outcome data
dat <- harmonise_data(exposure_dat, outcome_dat)

# Perform MR
mr_ast_sec <- mr(dat)

### MR effect of SEC31B expression on AST ###
# Get instruments from eqtlgen
exposure_dat <- extract_instruments("eqtl-a-ENSG00000075826")

# Get effects of instruments on outcome
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ukb-d-30650_irnt")

# Harmonise the exposure and outcome data
dat <- harmonise_data(exposure_dat, outcome_dat)

# Perform MR
mr_sec_ast <- mr(dat)

### MR effect of SEC31B expression in liver on AST ###
# Get instruments from gtex
gtex <- fread("/mnt/storage/private/mrcieu/data/broad/public/gtex/released/2020-03-09/data/v8_eQTL_all_associations/Liver.allpairs.txt.gz")
gtex <- gtex[startsWith(gtex$gene_id, "ENSG00000075826."),]
gtex <- gtex[gtex$pval_nominal < 0.01 & gtex$maf >= 0.05 & gtex$maf <= 0.95]
anno <- fread("https://storage.googleapis.com/gtex_analysis_v8/reference/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz")
gtex <- merge(gtex, anno, "variant_id")
exposure_dat <- format_data(
    gtex, 
    type="exposure", 
    eaf_col = "maf",
    snp_col = "rs_id_dbSNP151_GRCh38p7",
    beta_col = "slope",
    se_col = "slope_se",
    effect_allele_col = "alt",
    other_allele_col = "ref",
    phenotype_col= "gene_id",
    pval_col = "pval_nominal"
)
exposure_dat <- clump_data(exposure_dat)

# Get effects of instruments on outcome
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ukb-d-30650_irnt")

# Harmonise the exposure and outcome data
dat <- harmonise_data(exposure_dat, outcome_dat)

# Perform MR
mr_sec_gtex_ast <- mr(dat)

### LoF effect of SEC31B expression on AST ###

# read in extracted phenotypes
pheno <- fread("data/aspartate_aminotransferase.30650.txt")

# read in LoF calls
calls <- fread("data/ukbb.SEC31B.exome.txt.raw")

# link to IEU
linker <- get_filtered_linker(drop_standard_excl=TRUE, drop_non_white_british=TRUE, drop_related=TRUE, app="15825")
calls <- merge(calls, linker, by.x="FID", by.y="app15825")

# load variant annotations
anno <- fread("data/ukbb.SEC31B.exome.vep.txt", skip=42, check=T)
anno$stop <- grepl("frameshift_variant|stop_gained", anno$Consequence)

# select variants
high_impact <- anno %>%
    filter(Feature == "ENST00000370345") %>%
    filter(IMPACT == "HIGH") %>%
    select("X.Uploaded_variation") %>%
    distinct() %>%
    pull()

stop <- anno %>%
    filter(Feature == "ENST00000370345") %>%
    filter(stop) %>%
    select("X.Uploaded_variation") %>%
    distinct() %>%
    pull()

high_impact_missense <- anno %>%
    filter(Feature == "ENST00000370345") %>%
    filter(IMPACT == "HIGH" | (grepl("missense", Consequence) & SIFT == "deleterious" & PolyPhen == "probably_damaging")) %>%
    select("X.Uploaded_variation") %>%
    distinct() %>%
    pull()

# select high impact calls from total callset
high_impact_calls <- calls[,grepl(paste(high_impact, collapse="|"), names(calls)), with=F]
high_impact_calls <- data.frame(appieu=calls$appieu, n_high_impact=rowSums(high_impact_calls, na.rm=T))
high_impact_calls$has_high_impact <- high_impact_calls$n_high_impact != 0

stop_calls <- calls[,grepl(paste(stop, collapse="|"), names(calls)), with=F]
stop_calls <- data.frame(appieu=calls$appieu, n_stop=rowSums(stop_calls, na.rm=T))
stop_calls$has_stop <- stop_calls$n_stop != 0

high_impact_missense_calls <- calls[,grepl(paste(high_impact_missense, collapse="|"), names(calls)), with=F]
high_impact_missense_calls <- data.frame(appieu=calls$appieu, n_high_impact_missense=rowSums(high_impact_missense_calls, na.rm=T))
high_impact_missense_calls$has_high_impact_missense <- high_impact_missense_calls$n_high_impact_missense != 0

# merge
pheno <- merge(pheno, high_impact_calls, "appieu")
pheno <- merge(pheno, high_impact_missense_calls, "appieu")
pheno <- merge(pheno, stop_calls, "appieu")

# LoF effect on AST mean
fits <- lm(aspartate_aminotransferase.30650 ~ has_stop + age_at_recruitment.21022.0.0 + sex.31.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=pheno)
fitm <- lm(aspartate_aminotransferase.30650 ~ has_high_impact_missense + age_at_recruitment.21022.0.0 + sex.31.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=pheno)
fith <- lm(aspartate_aminotransferase.30650 ~ has_high_impact + age_at_recruitment.21022.0.0 + sex.31.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=pheno)