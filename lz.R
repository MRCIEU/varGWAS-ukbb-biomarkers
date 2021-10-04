library("data.table")
library("ieugwasr")
library("stringr")
library("dplyr")
source("funs.R")
set.seed(234)

# NOTE assumes LocusZoom and PLINK are on PATH

locuszoom <- function(gchr, gpos, ga1, ga2, gpval, chr, start, end, trait, refsnp){
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
  system(paste0("locuszoom --refsnp ", refsnp, " --cache None --rundir ./data --metal ", t ," -p ", trait, " --plotonly --pop EUR --build hg19 --source 1000G_Nov2014 --no-date --chr ", chr, " --start " ,start, " --end ", end))
}
locuszoom_gxe <- function(f, trait, label_file, chr, start, end, refsnp, main=T){
  d <- fread(f)
  d$p.value[which(d$p.value == 0)] <- .Machine$double.xmin
  t <- tempfile()
  d %>% 
    filter(main != grepl(":", term)) %>% 
    filter(grepl("^rs", term)) %>% 
    rowwise() %>% 
    mutate(MarkerName=str_split(term, "_", simplify=T)[1]) %>% 
    select(MarkerName, p.value) %>% 
    rename('P-value'="p.value") %>% 
    write.table(., sep="\t", row.names=F, quote=F, file=t)
  system(paste0("locuszoom --denote-markers-file ", label_file, " --refsnp ", refsnp, " --cache None --rundir ./data --metal ", t ," -p ", trait, " --plotonly --pop EUR --build hg19 --source 1000G_Nov2014 --no-date --chr ", chr, " --start " ,start, " --end ", end))
  snps <- fread(t, col.names=c("rsid", "pval"))
  snps <- ieugwasr::ld_clump(snps %>% filter(pval < 5e-8))
  print(snps)
}

# Urate stratified by sex
males <- fread("data/30880_irnt.gwas.imputed_v3.male.varorder.tsv")
males <- cbind(males, stringr::str_split(males$variant, ":", simplify=T) %>% as.data.frame(., stringsAsFactors=F) %>% dplyr::rename(chr="V1", pos="V2", ref="V3", alt="V4"))
females <- fread("data/30880_irnt.gwas.imputed_v3.female.varorder.tsv")
females <- cbind(females, stringr::str_split(females$variant, ":", simplify=T) %>% as.data.frame(., stringsAsFactors=F) %>% dplyr::rename(chr="V1", pos="V2", ref="V3", alt="V4"))
vgwas <- get_variants("urate.30880.0.0")
locuszoom(males$chr, males$pos, males$ref, males$alt, males$pval, "4", 10402838 - 500000, 10402838 + 500000, "Urate_M", "rs4530622")
locuszoom(females$chr, females$pos, females$ref, females$alt, females$pval, "4", 10402838 - 500000, 10402838 + 500000, "Urate_F", "rs4530622")
locuszoom(vgwas$chr, vgwas$pos, vgwas$oa, vgwas$ea, vgwas$phi_p, "4", 10402838 - 500000, 10402838 + 500000, "Urate_vGWAS", "rs4530622")

# GxE 1Mb around vQTL
locuszoom_gxe(
  "data/hdl_cholesterol.30760.0.0.x.sex.31.0.0.rs75627662.txt",
  "main_hdl_cholesterol.30760.0.0",
  "data/hdl_cholesterol.30760.0.0.x.sex.31.0.0.label.txt",
  "19", 45413576 - 500000, 45413576 + 500000,
  "rs75627662",
  main=T
)
locuszoom_gxe(
  "data/hdl_cholesterol.30760.0.0.x.sex.31.0.0.rs75627662.txt",
  "x_sex.31.0.0_hdl_cholesterol.30760.0.0",
  "data/hdl_cholesterol.30760.0.0.x.sex.31.0.0.label.txt",
  "19", 45413576 - 500000, 45413576 + 500000,
  "rs75627662",
  main=F
)

locuszoom_gxe(
  "data/triglycerides.30870.0.0.x.body_mass_index.21001.0.0.rs738409.txt",
  "main_triglycerides.30870.0.0",
  "data/triglycerides.30870.0.0.x.body_mass_index.21001.0.0.label.txt",
  "22", 44324727 - 500000, 44324727 + 500000,
  "rs738409",
  main=T
)
locuszoom_gxe(
  "data/triglycerides.30870.0.0.x.body_mass_index.21001.0.0.rs738409.txt",
  "x_body_mass_index.21001.0.0_triglycerides.30870.0.0",
  "data/triglycerides.30870.0.0.x.body_mass_index.21001.0.0.label.txt",
  "22", 44324727 - 500000, 44324727 + 500000,
  "rs738409",
  main=F
)

locuszoom_gxe(
  "data/urate.30880.0.0.x.sex.31.0.0.rs4530622.txt",
  "main_urate.30880.0.0",
  "data/urate.30880.0.0.x.sex.31.0.0.label.txt",
  "4", 10402838 - 500000, 10402838 + 500000,
  "rs4530622",
  main=T
)
locuszoom_gxe(
  "data/urate.30880.0.0.x.sex.31.0.0.rs4530622.txt",
  "x_sex.31.0.0_urate.30880.0.0",
  "data/urate.30880.0.0.x.sex.31.0.0.label.txt",
  "4", 10402838 - 500000, 10402838 + 500000,
  "rs4530622",
  main=F
)