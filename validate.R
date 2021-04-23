library('data.table')
library('optparse')
library('dplyr')
library('ieugwasr')
library('quantreg')
library('broom')
library("robustbase")
source("funs.R")
set.seed(124)

option_list = list(
  make_option(c("-p", "--pheno_file"), type="character", default=NULL, help="UKBB pheno CSV", metavar="character"),
  make_option(c("-t", "--trait"), type="character", default=NULL, help="Name of trait", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, help="Output file", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

model <- function(data, out, chr, pos, oa, ea, rsid) {
  # prepare data
  dosage <- extract_variant_from_bgen(chr, pos, oa, ea)
  data <- merge(data, dosage, "appieu")
  data <- na.omit(data)
  s <- paste0("chr", chr, "_", pos, "_", oa, "_", ea)
  s <- gsub(" ", "", s, fixed = TRUE)
  data$xsq <- (data %>% pull(!!s))^2

  # first-stage model
  f <- paste0(out, " ~ ", s, " + sex.31.0.0 + age_at_recruitment.21022.0.0 +", paste0("PC", seq(1, 10), collapse="+"))
  fit1q <- rq(as.formula(f), tau=0.5, data=data) # quantile
  fit1r <- tidy(lmrob(as.formula(f), data=data)) # SE robust

  # second-stage model
  data$d <- resid(fit1q)^2
  f <- paste0("d ~ ", s, " + xsq")
  fit2 <- lm(as.formula(f), data=data)

  # F-test
  fit0 <- lm(d ~ 1, data=data)
  ftest <- tidy(anova(fit0, fit2))
  fit2t <- tidy(fit2)

  return(data.frame(
      rsid=rsid,
      SNP=s,
      BETA_x=fit2t$estimate[2],
      BETA_xq=fit2t$estimate[3], 
      SE_x=fit2t$std.error[2],
      SE_xq=fit2t$std.error[3],
      Pvar=ftest$p.value[2],
      BETA=fit1r$estimate[2],
      SE=fit1r$std.error[2],
      Pmu=fit1r$p.value[2]
    )
  )
}

# read in extracted phenotypes
pheno <- fread(opt$p)

# load vGWAS for biomarker risk factor
gwas <- data.frame()
for (chr in seq(1,22)){
    if (chr < 10){
        file <- paste0("data/", opt$trait, ".vgwas.chr0", chr, ".txt")
    } else {
        file <- paste0("data/", opt$trait, ".vgwas.chr", chr, ".txt")
    }
    gwas <- rbind(gwas, fread(file))
}

# drop HLA region
gwas <- gwas[!(gwas$CHR == 6 & gwas$POS >= 29691116 & gwas$POS <= 33054976),]

# drop failed rows
gwas <- gwas %>%
    filter(SE_x != -1)

# select tophits
sig <- gwas[gwas$P < (5e-8/30) & gwas$EAF >= 0.05 & gwas$EAF <= 0.95]
sig <- sig[,c("RSID", "P")]
names(sig) <- c("rsid", "pval")
sig <- ld_clump(sig)

# subset rows
sig <- gwas[gwas$RSID %in% sig$rsid]

# GWAS
results <- apply(sig, 1, function(snp) {
  model(pheno, opt$trait, as.character(snp[['CHR']]), as.numeric(snp[['POS']]), as.character(snp[['OA']]), as.character(snp[['EA']]), as.character(snp[['RSID']]))
})
results <- rbindlist(results[!is.na(results)], fill=T)

# save assoc
write.csv(results, file=opt$o, row.names=F)