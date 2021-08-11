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

### MR of gene expression on biomarker ###

results <- results.interval %>% select(trait.vqtl, id.ensembl) %>% rename(trait.vqtl="trait")
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

# keep results with MR evidence
results.mr <- results.mr %>% filter(method %in% c("Inverse variance weighted", "Wald ratio")) %>% filter(pval < (0.05/n()))

# append IDs
results.mr$ensembl <- str_split(results.mr$id.exposure, "-", simplify=T)[,3]
results.mr$trait <- gsub("ukb-d-", "", results.mr$id.outcome)
results.mr$trait_id <- gsub("_irnt", "", results.mr$trait)
for (i in 1:nrow(results.mr)){
    results.mr$trait[i] <- grep(results.mr$trait_id[i], unique(results$trait), value = T)
}
results.mr$pair <- paste0(results.mr$ensembl, "-", results.mr$trait)

# subset interval data with MR evidence
results.interval$pair <- paste0(results.interval$id.ensembl, "-", results.interval$trait.vqtl)
results.interval.mr <- merge(results.interval, results.mr, "pair")

# unique
results.interval.mr$molecule_name.chembl <- NULL
results.interval.mr <- unique(results.interval.mr)

# write to table
results.interval.mr$effect <- paste0(round(results.interval.mr$b, 2), " (CI ", round(results.interval.mr$b - (1.96 * results.interval.mr$se), 2), ", ", round(results.interval.mr$b + (1.96 * results.interval.mr$se), 2), ")")
results.interval.mr %>% select(component_synonym.chembl, trait.vqtl, effect, rsid.vqtl, phi_p.vqtl) %>% unique