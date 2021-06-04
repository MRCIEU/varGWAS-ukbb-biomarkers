library("data.table")
library("dplyr")
library("DBI")
library("RSQLite")
library("multcomp")
library("ieugwasr")
set.seed(123)

# load data
pqtl <- fread("data/pqtls.csv")
pqtl <- pqtl %>% filter(eaf > 0.05 & eaf < 0.95)
mvqtl <- fread("data/vgwas.pqtl.txt")

# filter out non-mean effects
mvqtl <- mvqtl %>% filter(p < 0.1)

# merge data
dat <- merge(pqtl, mvqtl, "rsid", suffixes=c(".pqtl", ".mvqtl"))

# filter vQTLs
top <- dat[dat$phi_p < (5e-8/30)]

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

# filter mvQTLs at drug target loci
drugs <- merge(top, dat, by.x="gene", by.y="component_synonym")

# TODO compare RCT variance estimate with SNP estimate

### PCSK9 effect on LDL-c ###
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4845239/
# RCT - PCSK9 inhibitors lower LDL-c mean and variance (p=0.005; n=29)
# MR - rs191448950:A at PCSK9 lowers LDL-c mean and variance (p=3.15 x 10-8)