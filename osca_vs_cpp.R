# compare tophits with OSCA tools performed by Wang et al. Genotype-by-environment interactions inferred from genetic effects on phenotypic variability in the UK Biobank

# load summary stats
osca <- fread("BMI.ma")
gwas <- data.frame()
for (chr in seq(1,22)){
    if (chr < 10){
        path <- paste0("../jlst-cpp-vgwas/data//ukb_bmi.vgwas.chr0", chr, ".txt")
    } else {
        path <- paste0("../jlst-cpp-vgwas/data//ukb_bmi.vgwas.chr", chr, ".txt")
    }
    gwas <- rbind(gwas, fread(path))
}

# merge datasets
m < -merge(osca[,c("SNP", "p"),], gwas[,c("RSID", "P"),], by.x="SNP", by.y="RSID", suffixes=c(".osca", ".cpp"))

# correlate P values
cor.test(m$p, m$P)

#
#	Pearson's product-moment correlation
#
#data:  m$p and m$P
#t = 730.56, df = 5554547, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.2953215 0.2968390
#sample estimates:
#      cor 
#0.2960804 

# count n genome-wide hits in CPP where OSCA is genome-wide
table(m[m$p < 5e-8]$P < 5e-8)

#FALSE  TRUE 
# 2247  1947 

# count n genome-wide hits in OSCA where CPP is genome-wide
table(m[m$P < 5e-8]$p < 5e-8)

#FALSE  TRUE 
#  190  1947 

# count total number of genome-wide hits in OSCA
table(m$p < 5e-8)

#  FALSE    TRUE 
#5550355    4194 

# count total number of genome-wide hits in CPP
table(m$P < 5e-8)

#  FALSE    TRUE 
#5552412    2137