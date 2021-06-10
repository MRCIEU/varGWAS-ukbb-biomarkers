library("data.table")
library("dplyr")
library("stringr")
set.seed(1234)

# load gxg hits
snps <- fread("data/aspartate_aminotransferase.30650.0.0.gxg.txt")

# filter on p value
snps <- snps %>% filter(p.value < 5e-8)

# split term and make unique list of variants
snps <- cbind(snps, str_split(snps$term, ":", simplify=T), stringsAsFactors=F)
usnps <- unique(c(snps$V1,snps$V2), stringsAsFactors=F)
usnps <- as.data.frame(str_split(usnps, "_", simplify=T), stringsAsFactors=F)
usnps$V1 <- gsub("chr", "", usnps$V1)
usnps$V2 <- as.numeric(usnps$V2)

# write table to file
write.table(usnps, quote=F, row.names=F, sep="\t", file="data/aspartate_aminotransferase.30650.0.0.gxg-exome.txt")

### ON BC3 ###
gxg <- fread("aspartate_aminotransferase.30650.0.0.gxg-exome.txt")
gxg$key1 <- paste0(gxg$chr, ":", gxg$pos, ":", gxg$oa, ":", gxg$ea)
gxg$key2 <- paste0(gxg$chr, ":", gxg$pos, ":", gxg$ea, ":", gxg$oa)

results <- data.frame()
for (i in seq(1, 22)){
    bim <- fread(paste0("/projects/MRC-IEU/research/data/ukbiobank/genetic/variants/exome/dev/release_candidate/data/raw_downloaded/ukb_snp_chr",i,".bim"))
    bim <- bim %>% filter(V2 %in% gxg$key1 || V2 %in% gxg$key2)
    if (nrow(bim) > 1){
        results <- rbind(results, bim)
    }
}