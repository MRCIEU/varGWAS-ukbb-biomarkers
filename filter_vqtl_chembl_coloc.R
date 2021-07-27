library("DBI")
library("RSQLite")
library("dplyr")
library("data.table")
library("GenomicRanges")
library("IRanges")
library("stringr")
library("TwoSampleMR")
library("ieugwasr")
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

# add genomic coordinates for gene targets to chembl
ens <- fread("data/Homo_sapiens.GRCh37.82.bed")
names(ens)<-c("chrom", "start", "end", "ensg", "symbol")
anno <- merge(dat, ens, by.y="symbol", by.x="component_synonym")

# pad drug target loci by 500kb
anno$interval <- paste0(anno$chrom, ":", anno$start - 500000, "-", anno$end + 500000)
chembl_intervals <- GRanges(anno$interval)
stopifnot(length(chembl_intervals)==length(anno$interval))

### intersect colocalized regions drug target loci ###

# read coloc data
coloc <- fread("data/coloc.txt")

# select high prob of shared causal variant
coloc <- coloc[coloc$PP.H4.abf > .8]

# keep only eQTLgen data
coloc <- coloc %>% filter(grepl("^eqtl-a-", coloc$gene))

# split out gene from opengwas ID
coloc$ensg.coloc <- str_split(coloc$gene, "-", simplify=T)[,3]

# loop over coloc results and intersect with drug targets
results <- data.frame()
for (i in 1:nrow(coloc)){
    coloc_interval <- GRanges(coloc$region[i])
    
    # count overlaps between mvQTL and drug target
    counts <- suppressWarnings(countOverlaps(chembl_intervals, coloc_interval, type="any"))

    # save res
    if (sum(counts) > 0){
        res <- anno[counts>0,]
        res <- cbind(res, coloc[i,], stringsAsFactors=F)
        results <- rbind(results, res)
    }
}