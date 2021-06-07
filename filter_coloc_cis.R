library("data.table")
library("TwoSampleMR")
library("dplyr")
library("IRanges")
library("stringr")
library("GenomicRanges")
set.seed(123)

### Filter coloc results for evidence of shared casual variant at gene coding cis region ###

# read coloc data
coloc <- fread("data/coloc.txt")

# select high prob of shared causal variant
coloc <- coloc[coloc$PP.H4.abf > .8]

# append opengwas trait name to coloc data
ao <- available_outcomes()
ao <- ao %>% select(c("trait", "id"))
coloc <- merge(ao, coloc, by.x="id", by.y="gene")

# append sun et al target info to coloc data
lookup <- fread("data/001_SOMALOGIC_GWAS_protein_info.csv")
lookup <- lookup %>% select("Target", "TargetFullName")
coloc <- merge(coloc, lookup, by.x="trait.x", by.y="TargetFullName", all.x=T)

# append gene name to ENSG targets
bed <- fread("data/Homo_sapiens.GRCh37.82.bed")
bed <- bed %>% select("V4", "V5")
coloc <- merge(coloc, bed, by.x="trait.x", by.y="V4", all.x=T)

# carry over eQTL gene names
coloc$Target <- apply(coloc,1,function(x) if( is.na(x[["Target"]]) ) x[["V5"]] else x[["Target"]]   )

# manually add in Folkersen et al targets
coloc[which(coloc$trait.x=="KIT ligand"),]$Target <- "KITLG"
coloc[which(coloc$trait.x=="chitinase 3 like 1"),]$Target <- "CHI3L1"
coloc[which(coloc$trait.x=="coagulation factor III, tissue factor"),]$Target <- "F3"
coloc[which(coloc$trait.x=="follistatin"),]$Target <- "FST"
coloc[which(coloc$trait.x=="galectin 3"),]$Target <- "LGALS3"
coloc[which(coloc$trait.x=="platelet and endothelial cell adhesion molecule 1"),]$Target <- "PECAM1"
coloc[which(coloc$trait.x=="selectin E"),]$Target <- "SELE"

# add coordinates to gene targets
bed <- fread("data/Homo_sapiens.GRCh37.82.bed")
names(bed) <- c("chrom","start","end","ensg","gene")
bed$gene_region <- paste0(bed$chrom, ":", bed$start, "-", bed$end)
bed <- bed %>% select("gene","gene_region")
coloc <- merge(coloc, bed, by.x="Target", by.y="gene")
coloc$V5 <- NULL

# filter out trans regions
results <- data.frame()
for (i in 1:nrow(coloc)){
    # coloc region
    cregion <- GRanges(coloc$region[i])

    # gene region
    gregion <- GRanges(coloc$gene_region[i])

    # count overlap
    counts <- suppressWarnings(countOverlaps(cregion, gregion))

    # save res
    if (sum(counts) > 0){
        results <- rbind(results, coloc[i,])
    }
}

# save results
write.table(results, file="data/filter_coloc_cis.txt", sep="\t", quote=F, row.names=F)