library("ieugwasr")
library("dplyr")
source("funs.R")
set.seed(1234)
options(ieugwasr_api="http://64.227.44.193:8006/")

# load instruments for lipids at drug target loci & cross ref with vGWAS

# LDL-c

# get instruments from non-UKBB
ldl <- tophits("ieu-a-300", clump=0)

# filter drug target loci
hmgcr <- ldl %>% filter(chr == "5" & position > 74632154-500000 & position < 74657929+500000 & eaf >= 0.1 & eaf <= 0.9) %>% select(rsid, p) %>% rename(pval=p) %>% ld_clump
hmgcr <- ldl %>% filter(rsid %in% hmgcr$rsid) %>% mutate(gene="HMGCR")
pcsk9 <- ldl %>% filter(chr == "1" & position > 55505221-500000 & position < 55530525+500000 & eaf >= 0.1 & eaf <= 0.9) %>% select(rsid, p) %>% rename(pval=p) %>% ld_clump
pcsk9 <- ldl %>% filter(rsid %in% pcsk9$rsid) %>% mutate(gene="PCSK9")
npc1l1 <- ldl %>% filter(chr == "7" & position > 44552134-500000 & position < 44580914+500000 & eaf >= 0.1 & eaf <= 0.9) %>% select(rsid, p) %>% rename(pval=p) %>% ld_clump
npc1l1 <- ldl %>% filter(rsid %in% npc1l1$rsid) %>% mutate(gene="NPC1L1")
cetp <- ldl %>% filter(chr == "16" & position > 56995762-500000 & position < 57017757+500000 & eaf >= 0.1 & eaf <= 0.9) %>% select(rsid, p) %>% rename(pval=p) %>% ld_clump
cetp <- ldl %>% filter(rsid %in% cetp$rsid) %>% mutate(gene="CETP")

# orient effects to be lipid lowering
ldl <- rbind(hmgcr, pcsk9, npc1l1, cetp)
ldl$flip <- ldl$beta > 0
ldl$allele <- ""
ldl[which(ldl$flip),]$beta <- ldl[which(ldl$flip),]$beta * -1
ldl[which(ldl$flip),]$eaf <- 1 - ldl[which(ldl$flip),]$eaf
ldl[which(ldl$flip),]$allele <- ldl[which(ldl$flip),]$nea
ldl[which(ldl$flip),]$nea <- ldl[which(ldl$flip),]$ea
ldl[which(ldl$flip),]$ea <- ldl[which(ldl$flip),]$allele

# test for SNP effect on LDL-c variance


# HDL-c

# get instruments
hdl <- tophits("ieu-a-299", clump=0)

# filter drug target loci
cetp <- hdl %>% filter(chr == "16" & position > 56995762-500000 & position < 57017757+500000 & eaf >= 0.1 & eaf <= 0.9) %>% select(rsid, p) %>% rename(pval=p) %>% ld_clump
cetp <- hdl %>% filter(rsid %in% cetp$rsid) %>% mutate(gene="CETP")

# orient effects to be lipid lowering
hdl <- cetp
hdl$flip <- hdl$beta < 0
hdl$allele <- ""
hdl[which(hdl$flip),]$beta <- hdl[which(hdl$flip),]$beta * -1
hdl[which(hdl$flip),]$eaf <- 1 - hdl[which(hdl$flip),]$eaf
hdl[which(hdl$flip),]$allele <- hdl[which(hdl$flip),]$nea
hdl[which(hdl$flip),]$nea <- hdl[which(hdl$flip),]$ea
hdl[which(hdl$flip),]$ea <- hdl[which(hdl$flip),]$allele