library('optparse')
library('data.table')
library('broom')
source("funs.R")
set.seed(1234)

option_list = list(
  make_option(c("-p", "--pheno_file"), type="character", default=NULL, help="UKBB pheno CSV", metavar="character"),
  make_option(c("-t", "--trait"), type="character", default=NULL, help="Name of trait", metavar="character"),
  make_option(c("-s", "--snp_file"), type="character", default=NULL, help="SNP file", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, help="Output file", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

opt <- data.frame(p="data/hdl_cholesterol.30760.0.0.txt", s="data/hdl_cholesterol.30760.0.0.clump.txt", stringsAsFactors=F)

# read in extracted phenotypes
pheno <- fread(opt$p)

# read in snp list
snps <- fread(opt$s)

# load dosages
for (i in 1:nrow(snps)){
    dosage <- extract_variant_from_bgen(as.character(snps$chr[i]), as.double(snps$pos[i]), snps$oa[i], snps$ea[i])
    pheno <- merge(pheno, dosage, "appieu")
}

grep(names(pheno), "^chr")

# GWAS
results <- apply(snps, 1, function(snp) {
  mod(pheno, opt$t, as.character(snp[['chromosome']]), as.numeric(snp[['position']]), as.character(snp[['first_allele']]), as.character(snp[['alternative_alleles']]))
})
results <- rbindlist(results[!is.na(results)])

# save assoc
write.csv(results, file=opt$o, row.names=F)



pheno
names(pheno)
na.omit(pheno)
pheno <0 na.omit(pheno)
pheno <- na.omit(pheno)
# read in extracted phenotypes
pheno <- fread(opt$p)
snps <- fread(opt$s)
# load dosages
for (i in 1:nrow(snps)){
    dosage <- extract_variant_from_bgen(as.character(snps$chr[i]), as.double(snps$pos[i]), snps$oa[i], snps$ea[i])
    pheno <- merge(pheno, dosage, "appieu")
}
pheno
names(pheno)
grep(names(pheno), "^chr")
grep("^chr", names(pheno))
names(pheno)
names(pheno)
grep("^chr", names(pheno), value=T)
v <- grep("^chr", names(pheno), value=T)
grep("^chr", names(pheno), value=T)
grep("^chr", names(pheno), value=T) %>% filter(. != "chr16_56989590_T_C")
grep("^chr", names(pheno), value=T) %>% filter(., != "chr16_56989590_T_C")
grep("^chr", names(pheno), value=T) %>% filter("chr16_56989590_T_C")
filter(grep("^chr", names(pheno), value=T) %in% c("chr16_56989590_T_C"))
filter(grep("^chr", names(pheno), value=T) v %in% c("chr16_56989590_T_C"))
grep("^chr", names(pheno), value=T)
v <- grep("^chr", names(pheno), value=T)
v!= "chr16_56989590_T_C"
v[v!= "chr16_56989590_T_C"]
v <- v[v!= "chr16_56989590_T_C"]
paste0("chr16_56989590_T_C*",v)
paste0("chr16_56989590_T_C*",v, collapse="+)
paste0("chr16_56989590_T_C*",v, collapse="+")
paste0(paste0("chr16_56989590_T_C*",v, collapse="+"))
names(pheno)
paste0("age_at_recruitment.21022.0.0", paste0("chr16_56989590_T_C*",v, collapse="+"))
paste0("age_at_recruitment.21022.0.0+sex.31.0.0+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10", paste0("chr16_56989590_T_C*",v, collapse="+"))
paste0("hdl_cholesterol.30760.0.0 ~ age_at_recruitment.21022.0.0+sex.31.0.0+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+", paste0("chr16_56989590_T_C*",v, collapse="+"))
f<-as.formula(paste0("hdl_cholesterol.30760.0.0 ~ age_at_recruitment.21022.0.0+sex.31.0.0+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+", paste0("chr16_56989590_T_C*",v, collapse="+")))
f
lm(f, data=pheno)
summary(lm(f, data=pheno))
fit1 <- lm(f, data=pheno)
d2 <- resid(fit1)^2
lm(d2~ pheno$chr16_56989590_T_C)
lm(d2~ pheno$chr16_56989590_T_C)
f
f<-as.formula(paste0("d2 ~ age_at_recruitment.21022.0.0+sex.31.0.0+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+", paste0("chr16_56989590_T_C*",v, collapse="+")))
lm(f, data=cbind(pheno, d2))
summary(lm(f, data=cbind(pheno, d2)))
lm(d2~ pheno$chr16_56989590_T_C)
summary(lm(d2~ pheno$chr16_56989590_T_C))
summary(lm(f, data=cbind(pheno, d2)))
f<-as.formula(paste0("hdl_cholesterol.30760.0.0 ~ age_at_recruitment.21022.0.0+sex.31.0.0+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+", paste0("chr16_56989590_T_C*",v, collapse="+")))
lm(f, pheno)
summary(lm(f, pheno))
f
history()
history(max.show=9999)