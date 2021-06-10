library("data.table")
library("ieugwasr")
library("optparse")
library("dplyr")
library("coloc")
library("IRanges")
library("stringr")
library("GenomicRanges")
source("funs.R")
set.seed(123)

#H0: neither trait has a genetic association in the region
#H1: only trait 1 has a genetic association in the region.
#H2: only trait 2 has a genetic association in the region.
#H3: both traits are associated, but with different causal variants
#H4: both traits are associated, and share a single causal variant

# NOTE assumes LocusZoom and PLINK are on PATH

locuszoom <- function(gchr, gpos, ga1, ga2, gpval, chr, start, end, trait){
  message(paste0("Plotting region ", chr, ":", start, "-", end, " for trait ", trait))
  # subset data for interval
  gwas <- data.frame(gchr, gpos, ga1, ga2, gpval, stringsAsFactors=F)
  gwas <- gwas %>% filter(gchr==!!chr & gpos >=!!start & gpos <=!!end & nchar(ga1)==1 & nchar(ga2) ==1)
  gwas$MarkerName <- paste0("chr", gwas$gchr, ":", gwas$gpos)
  gwas <- gwas[,c("MarkerName", "gpval")]
  names(gwas)[2] <- "P-value"

  # write out records
  t <- tempfile()
  write.table(gwas, sep="\t", row.names=F, quote=F, file=t)

  # call LocusZoom
  system(paste0("locuszoom --cache None --rundir ./data --metal ", t ," -p ", trait, " --plotonly --pop EUR --build hg19 --source 1000G_Nov2014 --no-date --chr ", chr, " --start " ,start, " --end ", end))
}

option_list = list(
  make_option(c("-t", "--trait"), type="character", default=NULL, help="Name of trait", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, help="Output filename", metavar="character"),
  make_option(c("-s", "--snps"), type="character", default=NULL, help="File snp list", metavar="character"),
  make_option(c("-l", "--lz"), action="store_true", default=FALSE, help="Run locuszoom")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

message(paste0("working on ", opt$trait))

# load snps
mvqtl <- fread(opt$snps)
stopifnot(nrow(mvqtl)>0)

# load all vGWAS for trait
vgwas <- get_variants(opt$trait)

# find gene expression changes in blood associated with mvQTLs
genes <- phewas(mvqtl$rsid, pval = 5e-8, batch = c("eqtl-a", "prot-a", "prot-b", "prot-c"))
stopifnot(nrow(genes)>0)
ugenes <- unique(genes$id)

# drop mvQTLs that are not eQTLs
mvqtl <- mvqtl[mvqtl$rsid %in% genes$rsid]

# loop over mvQTLs for coloc
results <- data.frame()
for (i in 1:nrow(mvqtl)){
  message(paste0("Working on SNP: ", mvqtl$rsid[i]))

  # region to colocalise
  chr <- mvqtl$chr[i]
  pos <- mvqtl$pos[i]
  start <- pos - 500000
  end <- pos + 500000
  region <- paste0(chr, ":", start, "-", end)

  # subset vQTL data for interval
  vqtl <- vgwas[(vgwas$chr == chr & vgwas$pos >= start & vgwas$pos <= end)]

  for (j in 1:length(ugenes)){
    message(paste0("Working on expression of: ", ugenes[j]))

    # load eQTL data for interval
    eqtl <- associations(region, ugenes[j], proxies=0)

    if (nrow(eqtl) == 0){
      message(paste0("Skipping, not enough SNPs in interval"))
      next
    }

    if (!"eaf" %in% names(eqtl)){
      message(paste0("Skipping, eaf not given"))
      next
    }

    # exclude snps with missing values 
    eqtl <- eqtl[!is.na(eqtl$eaf),]
    eqtl <- eqtl[eqtl$eaf > 0,]
    eqtl <- eqtl[eqtl$eaf < 1,]

    if (nrow(eqtl) < 500){
      message(paste0("Skipping, not enough SNPs in interval"))
      next
    } else {
      message(paste0("Found N SNPs: ", nrow(eqtl)))
    }

    #format data for coloc
    gwas1.df <- list(
      N=vqtl$n,
      snp=vqtl$rsid,
      MAF=vqtl$eaf,
      pvalues=vqtl$phi_p,
      sdY=1,
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
    res$gene <- ugenes[j]
    res$region <- region
    res$trait <- opt$trait
    results <- rbind(results, res)
  }

}

# write out results
write.table(file=opt$o, results, sep="\t", quote=F, row.names=F)

if (opt$l){
  # select loci with shared casual variant
  h4 <- results[results$PP.H4.abf > .8,]
  df <- as.data.frame(str_split(h4$region, ":", 2, simplify=T), stringsAsFactors=F)
  names(df) <- c("chr", "interval")
  h4 <- cbind(h4, df)

  # collapse intervals for each gene and plot region
  for (gene in unique(h4$gene)){
    intervals <- h4[h4$gene == gene,]$interval

    # define overlapping interval
    intervals <- as.data.frame(
      GenomicRanges::reduce(
        GRanges(
          h4[h4$gene == gene,]$chr, 
          IRanges(
            h4[h4$gene == gene,]$interval
          )
        )
      )
    )

    # plot vGWAS
    locuszoom(vgwas$chr, vgwas$pos, vgwas$oa, vgwas$ea, vgwas$phi_p, intervals$seqnames, intervals$start, intervals$end, opt$trait)

    # extract eQTL
    eqtl <- associations(paste0(intervals$seqnames, ":", intervals$start, "-", intervals$end), gene, proxies=0)

    # plot eQTL
    locuszoom(eqtl$chr, eqtl$position, eqtl$ea, eqtl$nea, eqtl$p, intervals$seqnames, intervals$start, intervals$end, gene)
  }
}