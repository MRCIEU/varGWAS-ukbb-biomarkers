library("data.table")
library("TwoSampleMR")
library("ieugwasr")
library("dplyr")
library("IRanges")
library("stringr")
library("GenomicRanges")
source("funs.R")
set.seed(124)
options(ieugwasr_api="http://64.227.44.193:8006/")

# read coloc data
coloc <- fread("data/coloc.txt")

# select high prob of shared causal variant
coloc <- coloc[coloc$PP.H4.abf > .8]

# append opengwas trait name to coloc data
ao <- available_outcomes()
ao <- ao %>% select(c("trait", "id"))
coloc <- merge(coloc, ao, by.y="id", by.x="gene")

# append sun et al target info to coloc data
lookup <- fread("data/001_SOMALOGIC_GWAS_protein_info.csv")
lookup <- lookup %>% select("Target", "TargetFullName")
coloc <- merge(coloc, lookup, by.x="trait.y", by.y="TargetFullName", all.x=T)

# append gene name to ENSG targets
bed <- fread("data/Homo_sapiens.GRCh37.82.bed")
bed <- bed %>% select("V4", "V5")
coloc <- merge(coloc, bed, by.x="trait.y", by.y="V4", all.x=T)

# carry over eQTL gene names
coloc$Target <- apply(coloc,1,function(x) if( is.na(x[["Target"]]) ) x[["V5"]] else x[["Target"]]   )

# manually add in Folkersen et al targets
coloc[which(coloc$trait.y=="KIT ligand"),]$Target <- "KITLG"
coloc[which(coloc$trait.y=="chitinase 3 like 1"),]$Target <- "CHI3L1"
coloc[which(coloc$trait.y=="coagulation factor III, tissue factor"),]$Target <- "F3"
coloc[which(coloc$trait.y=="follistatin"),]$Target <- "FST"
coloc[which(coloc$trait.y=="galectin 3"),]$Target <- "LGALS3"
coloc[which(coloc$trait.y=="platelet and endothelial cell adhesion molecule 1"),]$Target <- "PECAM1"
coloc[which(coloc$trait.y=="selectin E"),]$Target <- "SELE"

# add neartest gene
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

# load vqtls
vqtl <- fread("data/vqtls.txt")
vqtl$key <- paste0(vqtl$chr, ":", vqtl$pos)

# merge
all <- merge(vqtl, ng, "key")
all$key <- paste0(all$chr, ":", all$pos, "-", all$pos)
all$coloc <- NA

# copy over
for (i in 1:nrow(all)){
    # snp region
    sregion <- GRanges(all$key[i])

    # gene region
    trait_coloc <- coloc %>% filter(trait.x == all$trait[i])
    gregion <- GRanges(trait_coloc %>% pull(region))

    # count overlap
    trait_coloc$counts <- suppressWarnings(countOverlaps(gregion, sregion))
    trait_coloc <- trait_coloc %>% filter(counts >0 ) %>% unique()

    if (nrow(trait_coloc) == 0){
        next
    }
    
    # add gene with coloc
    trait_coloc <- trait_coloc %>% select(Target, PP.H4.abf) %>% group_by(Target) %>% arrange(desc(PP.H4.abf)) %>% filter(row_number()==1)
    trait_coloc$gene <- paste0(trait_coloc$Target, " (H4=", round(trait_coloc$PP.H4.abf,2), ")")
    trait_coloc <- trait_coloc %>% arrange(desc(PP.H4.abf)) %>% head(n=3)
    f <- paste0(trait_coloc$gene, collapse=", ")

    # set field in maste table
    all$coloc[i] <- f
}

# combine gene col and only use nearest gene when no coloc evidence is available

# drop BMI
all <- all %>% filter(trait != "body_mass_index.21001.0.0")

# save table
all$Trait <- NA
for (i in 1:nrow(all)){
    all$Trait[i] <- biomarkers_abr[biomarkers==all$trait[i]]
}
all <- all %>% select(Trait, rsid, ea, gene, coloc, phi_p)
names(all) <- c("Trait", "RSID", "EA", "Nearest Gene", "Colocalized Gene", "P")

# order by trait and then by P
all <- all %>% arrange(Trait, P)

write.table(all, sep="\t", row.names=F, quote=F, file="all.vqtls.txt")