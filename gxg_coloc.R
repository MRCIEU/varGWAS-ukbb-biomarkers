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
usnps$v <- paste0(usnps$V1, ":", usnps$V2, "-", usnps$V2)

# find genes assocaited with gxg snps
genes <- phewas(usnps$v, pval = 5e-8, batch = c("eqtl-a", "prot-a", "prot-b", "prot-c"))
genes <- unique(genes)
stopifnot(nrow(genes)>0)

# pad regions for coloc
genes$region <- paste0(genes$chr, ":", genes$position - 250000, "-", genes$position + 250000)

results <- data.frame()
for (i in 1:nrow(genes)){
    # get assocs & drop duplicates
    trait <- associations(genes$region[i], opt$opengwas, proxies = 0)
    drop <- trait %>% select("rsid") %>% table %>% as.data.frame(., stringsAsFactors=F) %>% filter(Freq > 1) %>% pull(".")
    trait <- trait %>% filter(!rsid %in% drop)
    drop <- eqtl %>% select("rsid") %>% table %>% as.data.frame(., stringsAsFactors=F) %>% filter(Freq > 1) %>% pull(".")
    eqtl <- associations(genes$region[i], genes$id[i], proxies = 0)
    eqtl <- eqtl %>% filter(!rsid %in% drop)  

    #format data for coloc
    gwas1.df <- list(
      N=trait$n,
      beta=trait$beta,
      varbeta=trait$se^2,
      snp=trait$rsid,
      MAF=trait$eaf,
      type="quant"
    )

    gwas2.df <- list(
      N=eqtl$n,
      beta=eqtl$beta,
      varbeta=eqtl$se^2,
      snp=eqtl$rsid,
      MAF=eqtl$eaf,
      type="quant"
    )

    # perform coloc
    res <- coloc.abf(gwas1.df, gwas2.df)
    res <- as.data.frame(t(res$summary))
    res$gene <- genes$trait[i]
    res$position <- genes$position[i]
    res$chr <- genes$chr[i]
    res$trait <- opt$trait
    results <- rbind(results, res)
}

# select high evidence of shared casusl variant
results <- results %>% filter(PP.H4.abf > .7)
results$chrpos <- paste0(results$chr, ":", results$pos)

# merge on coloc data
results <- merge(results, snps, by.y="chrpos.1", by.x="chrpos", all.x=T)
results <- merge(results, snps, by.y="chrpos.2", by.x="chrpos", all.x=T)

# add gene to GxG data
snps <- snps %>% select(term, chrpos.1, chrpos.2, estimate, std.error, p.value)