library("DBI")
library("RSQLite")
library("dplyr")
library("data.table")
library("GenomicRanges")
library("IRanges")
library("stringr")
library("TwoSampleMR")
set.seed(123)

# extract tables
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

# merge
dat <- merge(MOLECULE_DICTIONARY, DRUG_MECHANISM, "molregno")
dat <- merge(dat, TARGET_DICTIONARY, "tid")
dat <- merge(dat, TARGET_COMPONENTS, "tid")
dat <- merge(dat, COMPONENT_SYNONYMS, "component_id", all.x=T)
dat <- merge(dat, COMPONENT_SEQUENCES, "component_id", all.x=T)
dat$component_id <- NULL
dat$tid <- NULL
dat$molregno <- NULL

# subset mvQTLs at drug target loci
d <- "/mnt/storage/scratch/ml18692/projects/jlst-cpp-vgwas/data"
files <- paste0(d, "/", list.files(d, pattern="*.clump.txt", full.names=F))
mvqtl <- data.frame()
for (file in files){
    f <- fread(file)
    f$trait <- gsub(".clump.txt", "", basename(file))
    mvqtl <- rbind(mvqtl, f)
}

# add genomic coordinates for gene targets
ens <- fread("data/Homo_sapiens.GRCh37.82.bed")
names(ens)<-c("chrom", "start", "end", "ensg", "symbol")
anno <- merge(dat, ens, by.y="symbol", by.x="component_synonym")
anno$interval <- paste0(anno$chrom, ":", anno$start - 500000, "-", anno$end + 500000)
chembl_intervals <- GRanges(anno$interval)
stopifnot(length(chembl_intervals)==length(anno$interval))

# intersect vQTLs at drug target loci
results <- data.frame()
for (i in 1:nrow(mvqtl)){
    mvqtl_interval <- GRanges(paste0(mvqtl$chr[i], ":", mvqtl$pos[i], "-", mvqtl$pos[i]))
    
    # count overlaps between mvQTL and drug target
    counts <- suppressWarnings(countOverlaps(chembl_intervals, mvqtl_interval, type="any"))

    # save res
    if (sum(counts) > 0){
        res <- anno[counts>0,]
        res <- cbind(res, mvqtl[i], stringsAsFactors=F)
        results <- rbind(results, res)
    }
}

# process eQTL
erels <- unique(results[,c("ensg", "trait")]) %>% filter(trait!="body_mass_index.21001.0.0")
names(erels)[1] <- "target"
erels$target <- paste0("eqtl-a-", erels$target)

# process pQTL
ao <- available_outcomes()
sun <- ao %>% filter(pmid == 29875488)
lookup <- fread("data/001_SOMALOGIC_GWAS_protein_info.csv")
sun <- merge(sun, lookup, by.x="trait", by.y="TargetFullName")
genes <- fread("data/soma_genes.txt", check.names=T)
genes <- genes %>% filter(Chromosome.scaffold.name %in% c("X", seq(1,22)))
sun <- merge(sun, genes, by.x="Target", by.y="Gene.name")
prels <- unique(results[,c("ensg", "trait")]) %>% filter(trait!="body_mass_index.21001.0.0")
prels <- merge(prels, sun, by.x="ensg", by.y="Gene.stable.ID") 
prels <- unique(prels[c("id", "trait.x")])
names(prels)[1] <- "target"
names(prels)[2] <- "trait"

# combine eQTL and pQTL
rels <- rbind(erels, prels)
rels$trait <- str_split(rels$trait, "\\.", simplify=T)[,2]
rels$trait <- paste0("ukb-d-", rels$trait, "_irnt")

# test for MR effect of eQTL/pQTL on trait mean & bi-directional effect
mr_results <- data.frame()
for (i in 1:nrow(rels)){
    # Test for effect of gene product on biomarker conc

    # Get instruments
    exposure_dat <- extract_instruments(rels$target[i])

    if (is.null(exposure_dat) || nrow(exposure_dat) == 0){
        next
    }

    # Get effects of instruments on outcome
    outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes=rels$trait[i])

    # Harmonise the exposure and outcome data
    dat <- harmonise_data(exposure_dat, outcome_dat)

    # Perform MR
    res <- mr(dat)

    # Steiger test
    result = tryCatch({
        directionality_test(dat)
    }, error = function(error_condition) {
        NULL
    })
    if (!is.null(result)){
      res <- cbind(res, result[5:8])
    }

    mr_results <- plyr::rbind.fill(mr_results, res)

    # Test for effect of biomarker on gene product

    # Get instruments
    exposure_dat <- extract_instruments(rels$trait[i])

    # Get effects of instruments on outcome
    outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes=rels$target[i])

    # Harmonise the exposure and outcome data
    dat <- harmonise_data(exposure_dat, outcome_dat)

    # Perform MR
    res <- mr(dat)

    # Steiger test
    result = tryCatch({
        directionality_test(dat)
    }, error = function(error_condition) {
        NULL
    })
    if (!is.null(result)){
      res <- cbind(res, result[5:8])
    }

    mr_results <- plyr::rbind.fill(mr_results, res)
}

# add key
mr_results$key <- paste0(mr_results$id.exposure, ".", mr_results$id.outcome)

# Steiger filter
mr_results <- mr_results[mr_results$correct_causal_direction,]

# drop null MR estimates for IVW / Wald
main <- mr_results %>% filter(method == "Inverse variance weighted" | method == "Wald ratio")
main <- main %>% filter(pval < (0.05/nrow(main)))

# drop estimates not supported by sensitivity analyses
sens <- mr_results %>% filter(method != "Inverse variance weighted" & method != "Wald ratio") %>% filter(pval < 0.05)
sig <- main[main$key %in% sens$key,]
sig <- rbind(sig, main[main$method == "Wald ratio",])

# filter out bidirectional effects
bi <- sig[grepl("irnt", sig$exposure),]
sig <- sig[!grepl("irnt", sig$exposure),]
if (nrow(bi) > 0){
    bi$rev_key <- paste0(bi$id.outcome, ".", bi$id.exposure)
    sig <- sig[!sig$key %in% bi$rev_key,]
}