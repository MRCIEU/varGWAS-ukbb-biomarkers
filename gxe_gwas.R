load("data/pheno.RData")
library('optparse')
library('data.table')
library('dplyr')
library('broom')
library('ieugwasr')
source("funs.R")
set.seed(1234)
options(ieugwasr_api="http://64.227.44.193:8006/")

option_list = list(
  make_option(c("-t", "--trait"), type="character", default=NULL, help="Name of trait", metavar="character"),
  make_option(c("-m", "--modifier"), type="character", default=NULL, help="Name of modifier", metavar="character"),
  make_option(c("-c", "--chr"), type="integer", default=NULL, help="Chromosome", metavar="character"),
  make_option(c("-s", "--start"), type="integer", default=NULL, help="Start", metavar="character"),
  make_option(c("-e", "--end"), type="integer", default=NULL, help="End", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Output", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# load linker
linker <- get_filtered_linker(drop_standard_excl=TRUE, drop_non_white_british=TRUE, drop_related=TRUE, application="15825")

# load covariates
covariates <- get_covariates()
pc <- get_genetic_principal_components()

# merge data
dat <- merge(linker, covariates, "appieu")
dat <- merge(dat, pheno, by.x="app15825", by.y="eid")
dat <- merge(dat, pc, "appieu")

# SD
dat[[opt$trait]] <- dat[[opt$trait]] / sd(dat[[opt$trait]], na.rm=T)
dat[[opt$modifier]] <- dat[[opt$modifier]] / sd(dat[[opt$modifier]], na.rm=T)

# load snps
if (opt$chr < 10){
    variants <- bgen.load(paste0("/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr0", opt$chr, ".bgen"),
        data.frame(chromosome = paste0("0", opt$chr), start = opt$start, end = opt$end)
    )
} else {
    variants <- bgen.load(paste0("/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr", opt$chr, ".bgen"),
        data.frame(chromosome = paste0(opt$chr), start = opt$start, end = opt$end)
    )
}
dosage <- as.data.frame(
    apply(variants$data, 1, function(data) { return(data[,1]*0 + data[,2]*1 + data[,3]*2) })
)
names(dosage) <- paste0(names(dosage), "_", variants$variants$allele0, "_", variants$variants$allele1)
snps <- names(dosage)
dosage$appieu <- row.names(dosage)

message(paste0("Found ", length(names(dosage))))

# merge geno & pheno
dat <- merge(dat, dosage, "appieu")

# test GxE
results <- data.frame()
for (snp in snps){
    f <- paste0(opt$trait, " ~ ", snp, " * ", opt$modifier, " + age_at_recruitment.21022.0.0 + sex.31.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")
    fit <- lm(as.formula(f), data=dat)
    fit <- tidy(fit)
    term <- fit %>% dplyr::filter(grepl(":", term))
    results <- rbind(results, term)
    term <- fit %>% dplyr::filter(term == !!snp)
    results <- rbind(results, term)
}

write.table(results, row.names=F, quote=F, sep="\t", file=opt$o)
