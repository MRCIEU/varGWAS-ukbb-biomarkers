library('optparse')
library('data.table')
library('broom')
set.seed(1234)

# Run B-P model using R to validate cpp application

option_list = list(
  make_option(c("-p", "--pheno_file"), type="character", default=NULL, help="UKBB pheno CSV", metavar="character"),
  make_option(c("-t", "--trait"), type="character", default=NULL, help="Name of trait", metavar="character"),
  make_option(c("-s", "--snp_file"), type="character", default=NULL, help="SNP file", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, help="Output file", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

mod <- function(snp, dat, model) {
    dosage <- extract_variant_from_bgen(snp[['chr']], as.numeric(snp[['pos']]), snp[['other_allele']], snp[['effect_allele']])
    fit <- lm(model, data=merge(dosage, dat, "appieu"))
    return(tidy(fit)[2,])
}

# read in extracted phenotypes
pheno <- fread(opt$p)

# read in snp list
snps <- fread(opt$s)

# GWAS
results <- apply(snps, 1, function(snp) {
    adj[1] <- paste0("chr", paste(snp[['chr']], snp[['pos']], snp[['other_allele']], snp[['effect_allele']], sep="_"))
    model <- paste(adj, collapse = " + ")
    lm_assoc(snp, dat, as.formula(paste0(opt$t, " ~ ", model)))
})
results <- rbindlist(results)

# save assoc
write.csv(results, file=opt$o, row.names=F)