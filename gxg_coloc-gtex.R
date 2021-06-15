library("data.table")
library('optparse')
library('dplyr')
library('ieugwasr')
library("stringr")
library('coloc')
source("funs.R")
set.seed(123)

option_list = list(
  make_option(c("-t", "--trait"), type="character", default=NULL, help="Name of trait", metavar="character"),
  make_option(c("-o", "--opengwas"), type="character", default=NULL, help="OpenGWAS ID", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# read in GxG associations
snps <- fread(paste0("data/", opt$trait, ".gxg.txt"))

# take top GxG hits for coloc
snps <- snps %>% filter(p.value < 0.05 / (1e+6 * 30 + 250000 * 30))
snps <- cbind(snps, str_split(snps$term, ":", simplify=T), stringsAsFactors=F)
v1 <- as.data.frame(str_split(snps$V1, "_", simplify=T), stringsAsFactors=F)
names(v1) <- c("chr.1", "pos.1", "oa.1", "ea.1")
v2 <- as.data.frame(str_split(snps$V2, "_", simplify=T), stringsAsFactors=F)
names(v2) <- c("chr.2", "pos.2", "oa.2", "ea.2")
snps <- cbind(snps, v1, v2, stringsAsFactors=F)
snps$chr.1 <- gsub("chr", "", snps$chr.1)
snps$chr.2 <- gsub("chr", "", snps$chr.2)
snps$pos.1 <- as.numeric(snps$pos.1)
snps$pos.2 <- as.numeric(snps$pos.2)
snps$chrpos.1 <- paste0(snps$chr.1, ":", snps$pos.1)
snps$chrpos.2 <- paste0(snps$chr.2, ":", snps$pos.2)

# split term
usnps <- unique(c(snps$V1,snps$V2), stringsAsFactors=F)
usnps <- as.data.frame(str_split(usnps, "_", simplify=T), stringsAsFactors=F)
usnps$V1 <- gsub("chr", "", usnps$V1)
usnps$V2 <- as.numeric(usnps$V2)

# load GTEx data
gtex <- fread("/mnt/storage/private/mrcieu/data/broad/public/gtex/released/2018-10-05/data/GTEx_Analysis_v7_eQTL_all_associations/Liver.allpairs.txt.gz")
gtex <- cbind(gtex, as.data.frame(str_split(gtex$variant_id, "_", simplify=T), stringsAsFactors=F))
gtex$V2 <- as.numeric(gtex$V2)

# find genes assocaited with gxg snps
results <- data.frame()
for (i in 1:nrow(usnps)){
  message(paste0("Working on: ", paste0(usnps$V1[i], "-", usnps$V2[i], "-", usnps$V3[i],"-", usnps$V4[i])))
  gtex_snps <- gtex %>% filter(V1 == usnps$V1[i] & V2 >= (usnps$V2[i] - 250000) & V2 <= (usnps$V2[i] + 250000))
  gtex_snps <- gtex_snps[complete.cases(gtex_snps)]
  if (gtex_snps %>% filter(pval_nominal < 5e-5) %>% nrow > 0){
    # get assocs & drop duplicates
    trait <- associations(paste0(usnps$V1[i], ":", usnps$V2[i] - 250000, "-", usnps$V2[i] + 250000), opt$opengwas, proxies = 0)
    if (min(trait$p) > 5e-5){
      next
    }
    drop <- trait %>% select("rsid") %>% table %>% as.data.frame(., stringsAsFactors=F) %>% filter(Freq > 1) %>% pull(".")
    trait <- trait %>% filter(!rsid %in% drop)
    trait$key <- paste0(trait$chr, "_", trait$position, "_", trait$nea, "_", trait$ea, "_b37")

    for (gene in unique(gtex_snps$gene_id)){
      message(paste0("Working on: ", gene))
      gtex_snps_tmp <- gtex_snps %>% filter(gene_id == !!gene)

      #format data for coloc
      gwas1.df <- list(
        N=trait$n,
        beta=trait$beta,
        varbeta=trait$se^2,
        snp=trait$key,
        MAF=trait$eaf,
        type="quant"
      )

      gwas2.df <- list(
        N=rep(200, nrow(gtex_snps_tmp)),
        beta=gtex_snps_tmp$slope,
        varbeta=gtex_snps_tmp$slope_se^2,
        snp=gtex_snps_tmp$variant_id,
        MAF=gtex_snps_tmp$maf,
        type="quant"
      )

      # perform coloc
      res <- coloc.abf(gwas1.df, gwas2.df)
      res <- as.data.frame(t(res$summary))
      res$gene <- gene
      res$position <- usnps$V2[i]
      res$chr <- usnps$V1[i]
      res$trait <- opt$trait
      results <- rbind(results, res)
    }

  }
}

# write out all gtex results
write.table(results, sep="\t", quote=F, row.names=F, file=paste0("data/", opt$trait, ".gxg-coloc-gtex.txt"))

# merge on GxG effects
