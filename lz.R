library("data.table")
library("dplyr")
source("funs.R")
set.seed(234)

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
  system(paste0("locuszoom --refsnp rs4530622 --cache None --rundir ./data --metal ", t ," -p ", trait, " --plotonly --pop EUR --build hg19 --source 1000G_Nov2014 --no-date --chr ", chr, " --start " ,start, " --end ", end))
}

males <- fread("data/30880_irnt.gwas.imputed_v3.male.varorder.tsv")
males <- cbind(males, stringr::str_split(males$variant, ":", simplify=T) %>% as.data.frame(., stringsAsFactors=F) %>% dplyr::rename(chr="V1", pos="V2", ref="V3", alt="V4"))
females <- fread("data/30880_irnt.gwas.imputed_v3.female.varorder.tsv")
females <- cbind(females, stringr::str_split(females$variant, ":", simplify=T) %>% as.data.frame(., stringsAsFactors=F) %>% dplyr::rename(chr="V1", pos="V2", ref="V3", alt="V4"))
vgwas <- get_variants("urate.30880.0.0")

locuszoom(males$chr, males$pos, males$ref, males$alt, males$pval, "4", 10402838 - 500000, 10402838 + 500000, "Urate_M")
locuszoom(females$chr, females$pos, females$ref, females$alt, females$pval, "4", 10402838 - 500000, 10402838 + 500000, "Urate_F")
locuszoom(vgwas$chr, vgwas$pos, vgwas$oa, vgwas$ea, vgwas$phi_p, "4", 10402838 - 500000, 10402838 + 500000, "Urate_vGWAS")