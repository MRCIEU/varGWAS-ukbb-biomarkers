library("data.table")
library("ieugwasr")
library("stringr")
library("dplyr")
source("funs.R")
set.seed(234)
options(ieugwasr_api="http://64.227.44.193:8006/")


# NOTE assumes LocusZoom and PLINK are on PATH

locuszoom <- function(gchr, gpos, ga1, ga2, gpval, chr, start, end, trait, refsnp){
  message(paste0("Plotting region ", chr, ":", start, "-", end, " for trait ", trait))
  # subset data for interval
  gwas <- data.frame(gchr, gpos, ga1, ga2, gpval, stringsAsFactors=F)
  gwas <- gwas %>% dplyr::filter(gchr==!!chr & gpos >=!!start & gpos <=!!end & nchar(ga1)==1 & nchar(ga2) ==1)
  gwas$MarkerName <- paste0("chr", gwas$gchr, ":", gwas$gpos)
  gwas <- gwas[,c("MarkerName", "gpval")]
  names(gwas)[2] <- "P-value"

  # write out records
  t <- tempfile()
  write.table(gwas, sep="\t", row.names=F, quote=F, file=t)

  # call LocusZoom
  system(paste0("locuszoom --refsnp ", refsnp, " --cache None --rundir ./data --metal ", t ," -p ", trait, " --plotonly --pop EUR --build hg19 --source 1000G_Nov2014 --no-date --chr ", chr, " --start " ,start, " --end ", end))
}
locuszoom_gxe <- function(f, trait, label_file, chr, start, end, refsnp, main=T){
  d <- fread(f)
  d$p.value[which(d$p.value == 0)] <- .Machine$double.xmin
  t <- tempfile()
  d %>% 
    dplyr::filter(main != grepl(":", term)) %>% 
    dplyr::filter(grepl("^rs", term)) %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(MarkerName=str_split(term, "_", simplify=T)[1]) %>% 
    dplyr::select(MarkerName, p.value) %>% 
    dplyr::rename('P-value'="p.value") %>% 
    write.table(., sep="\t", row.names=F, quote=F, file=t)
  system(paste0("locuszoom --denote-markers-file ", label_file, " --refsnp ", refsnp, " --cache None --rundir ./data --metal ", t ," -p ", trait, " --plotonly --pop EUR --build hg19 --source 1000G_Nov2014 --no-date --chr ", chr, " --start " ,start, " --end ", end))
  snps <- fread(t, col.names=c("rsid", "pval"))
  snps <- ieugwasr::ld_clump(snps %>% filter(pval < 5e-8))
  print(snps)
}

# Urate-SLC2A9-by-sex
males <- fread("data/30880_irnt.gwas.imputed_v3.male.varorder.tsv.gz")
males <- cbind(males, stringr::str_split(males$variant, ":", simplify=T) %>% as.data.frame(., stringsAsFactors=F) %>% dplyr::rename(chr="V1", pos="V2", ref="V3", alt="V4"))
males$pos <- as.numeric(males$pos)
females <- fread("data/30880_irnt.gwas.imputed_v3.female.varorder.tsv.gz")
females <- cbind(females, stringr::str_split(females$variant, ":", simplify=T) %>% as.data.frame(., stringsAsFactors=F) %>% dplyr::rename(chr="V1", pos="V2", ref="V3", alt="V4"))
females$pos <- as.numeric(females$pos)
vgwas <- get_variants("urate.30880.0.0")
locuszoom(males$chr, males$pos, males$ref, males$alt, males$pval, "4", 9920347 - 500000, 9920347 + 500000, "Urate_M", "rs10805346")
locuszoom(females$chr, females$pos, females$ref, females$alt, females$pval, "4", 9920347 - 500000, 9920347 + 500000, "Urate_F", "rs10805346")
locuszoom(vgwas$chr, vgwas$pos, vgwas$oa, vgwas$ea, vgwas$phi_p, "4", 9920347 - 500000, 9920347 + 500000, "Urate_var", "rs10805346")