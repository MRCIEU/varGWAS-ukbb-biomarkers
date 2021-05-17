# compare tophits with OSCA tools performed by Wang et al. Genotype-by-environment interactions inferred from genetic effects on phenotypic variability in the UK Biobank
library('data.table')
library('dplyr')
library('broom')
source("funs.R")
set.seed(1234)

# load summary stats
osca <- fread("../13-QT-vQTL-summary-data/BMI.ma")
gwas <- get_variants("body_mass_index.21001.0.0")

# merge datasets
m <- merge(osca[,c("SNP", "p"),], gwas[,c("rsid", "phi_p"),], by.x="SNP", by.y="rsid")

# correlate P values
tidy(cor.test(m$p, m$phi_p))

#  estimate statistic p.value parameter conf.low conf.high method     alternative
#     <dbl>     <dbl>   <dbl>     <int>    <dbl>     <dbl> <chr>      <chr>      
#1    0.378      960.       0   5512280    0.378     0.379 Pearson's~ two.sided  

# count number of OSCA hits that are also B-P hits
table(m[m$p < 5e-8]$phi_p < 5e-8)

#FALSE  TRUE 
# 1228  1295 

# count number of B-P hits that are also OSCA hits
table(m[m$phi_p < 5e-8]$p < 5e-8)

#FALSE  TRUE 
#   51  1295 

# count total number of genome-wide hits in OSCA
table(m$p < 5e-8)

#  FALSE    TRUE 
#5509759    2523 

# count total number of genome-wide hits in B-P
table(m$phi_p < 5e-8)

#  FALSE    TRUE 
#5510936    1346 

# is OSCA more sensitive than B-P or higher T1E?