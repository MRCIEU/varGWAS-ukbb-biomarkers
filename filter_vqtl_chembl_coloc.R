library("DBI")
library("RSQLite")
library("dplyr")
library("data.table")
library("GenomicRanges")
library("IRanges")
library("stringr")
library("TwoSampleMR")
library("ieugwasr")
library("stringr")
source("funs.R")
set.seed(124)
options(ieugwasr_api="http://64.227.44.193:8006/")

### prepare chembl data ###

# extract drug-target-gene data from chembl
con <- dbConnect(drv=RSQLite::SQLite(), dbname="data/chembl/chembl_28/chembl_28_sqlite/chembl_28.db")
MOLECULE_DICTIONARY <- dbGetQuery(conn=con, statement="SELECT molregno,pref_name,molecule_type FROM MOLECULE_DICTIONARY")
names(MOLECULE_DICTIONARY)[names(MOLECULE_DICTIONARY) == "pref_name"] <- "molecule_name"
DRUG_MECHANISM <- dbGetQuery(conn=con, statement="SELECT molregno,mechanism_of_action,tid,action_type FROM DRUG_MECHANISM")
TARGET_DICTIONARY <- dbGetQuery(conn=con, statement="SELECT tid,target_type,pref_name FROM TARGET_DICTIONARY WHERE organism=='Homo sapiens'")
names(TARGET_DICTIONARY)[names(TARGET_DICTIONARY) == "pref_name"] <- "target_name"
TARGET_COMPONENTS <- dbGetQuery(conn=con, statement="SELECT tid,component_id FROM TARGET_COMPONENTS")
COMPONENT_SYNONYMS <- dbGetQuery(conn=con, statement="SELECT component_id,component_synonym FROM COMPONENT_SYNONYMS WHERE syn_type=='GENE_SYMBOL'")
COMPONENT_SEQUENCES <- dbGetQuery(conn=con, statement="SELECT component_id,accession FROM COMPONENT_SEQUENCES WHERE organism=='Homo sapiens'")
dbDisconnect(con)

# merge records
dat <- merge(MOLECULE_DICTIONARY, DRUG_MECHANISM, "molregno")
dat <- merge(dat, TARGET_DICTIONARY, "tid")
dat <- merge(dat, TARGET_COMPONENTS, "tid")
dat <- merge(dat, COMPONENT_SYNONYMS, "component_id", all.x=T)
dat <- merge(dat, COMPONENT_SEQUENCES, "component_id", all.x=T)
dat$component_id <- NULL
dat$tid <- NULL
dat$molregno <- NULL
names(dat) <- paste0(names(dat), ".chembl")

# add genomic coordinates for gene targets to chembl
ens <- fread("data/Homo_sapiens.GRCh37.82.bed")
names(ens)<-c("chrom.ensembl", "start.ensembl", "end.ensembl", "id.ensembl", "symbol.ensembl")
anno <- merge(dat, ens, by.y="symbol.ensembl", by.x="component_synonym.chembl")

# pad drug target loci by 500kb
anno$interval.ensembl <- paste0(anno$chrom.ensembl, ":", anno$start.ensembl - 500000, "-", anno$end.ensembl + 500000)
chembl_intervals <- GRanges(anno$interval.ensembl)
stopifnot(length(chembl_intervals)==length(anno$interval.ensembl))

### intersect vQTL with drug target loci ###

vqtls <- fread("data/vqtls.txt") %>% select(trait, chr, pos, rsid, oa, ea, phi_p)
names(vqtls) <- paste0(names(vqtls), ".vqtl")
vqtls$interval.vqtl <- paste0(vqtls$chr.vqtl, ":", vqtls$pos.vqtl, "-", vqtls$pos.vqtl)

# loop over vQTLs and intersect with drug targets
results.interval <- data.frame()
for (i in 1:nrow(vqtls)){
    vqtl_interval <- GRanges(vqtls$interval.vqtl[i])
    
    # count overlaps between mvQTL and drug target
    counts <- suppressWarnings(countOverlaps(chembl_intervals, vqtl_interval, type="any"))

    # save res
    if (sum(counts) > 0){
        res <- anno[counts>0,]
        res <- cbind(res, vqtls[i,], stringsAsFactors=F)
        results.interval <- rbind(results.interval, res)
    }
}

### intersect colocalized regions drug target loci ###

# read coloc data
coloc <- fread("data/coloc.txt")
coloc <- coloc %>% filter(trait.coloc != "body_mass_index.21001.0.0")

# select high prob of shared causal variant
coloc <- coloc[coloc$PP.H4.abf > .8]

# keep only eQTLgen data
coloc <- coloc %>% filter(grepl("^eqtl-a-", coloc$gene))

# split out gene from opengwas ID
coloc$ensembl <- str_split(coloc$gene, "-", simplify=T)[,3]
coloc$gene <- NULL

# append coloc suffix
names(coloc) <- paste0(names(coloc), ".coloc")

# loop over coloc results and intersect with drug targets
results.coloc <- data.frame()
for (i in 1:nrow(coloc)){
    coloc_interval <- GRanges(coloc$region.coloc[i])
    
    # count overlaps between mvQTL and drug target
    counts <- suppressWarnings(countOverlaps(chembl_intervals, coloc_interval, type="any"))

    # save res
    if (sum(counts) > 0){
        res <- anno[counts>0,]
        res <- cbind(res, coloc[i,], stringsAsFactors=F)
        results.coloc <- rbind(results.coloc, res)
    }
}

# keep overlap of gene between coloc and ensembl
results.coloc <- results.coloc %>% filter(ensembl.coloc == id.ensembl)

### MR ###
results <- rbind(results.interval %>% select(trait.vqtl, id.ensembl) %>% rename(trait.vqtl="trait"), results.coloc %>% select(trait.coloc, id.ensembl) %>% rename(trait.coloc="trait"))
results <- unique(results)
results$id.outcome <- str_split(results$trait, "\\.", simplify=T)[,2] %>% paste0("ukb-d-", ., "_irnt")
results$id.exposure <- paste0("eqtl-a-", results$id.ensembl)

# drop BMI
results <- results %>% filter(id.outcome != "ukb-d-21001_irnt")

# test for MR evidence of gene expression on biomarker
results.mr <- data.frame()
for (i in 1:nrow(results)){
    exposure_dat <- extract_instruments(results$id.exposure[i])

    if (is.null(exposure_dat)){
        next
    }

    # Get effects of instruments on outcome
    outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes=results$id.outcome[i])

    # Harmonise the exposure and outcome data
    dat <- harmonise_data(exposure_dat, outcome_dat)

    # Perform MR
    res <- mr(dat)

    results.mr <- rbind(results.mr, res)
}