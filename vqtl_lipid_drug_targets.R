library("ieugwasr")
library("dplyr")
source("funs.R")
set.seed(1234)
options(ieugwasr_api="http://64.227.44.193:8006/")

# load instruments for lipids at drug target loci & cross ref with vGWAS

# LDL-c

# get instruments
ldl <- tophits("ieu-b-110", clump=0)

# load vGWAS for LDL-c
vgwas <- get_variants("ldl_direct.30780.0.0")

# filter drug target loci
hmgcr <- ldl %>% filter(chr == "5" & position > 74632154-500000 & position < 74657929+500000)
hmgcr <- vgwas %>% filter(rsid %in% hmgcr$rsid)
hmgcr$ gene <- "HMGCR"
pcsk9 <- ldl %>% filter(chr == "1" & position > 55505221-500000 & position < 55530525+500000)
pcsk9 <- vgwas %>% filter(rsid %in% pcsk9$rsid)
pcsk9$ gene <- "PCSK9"
npc1l1 <- ldl %>% filter(chr == "7" & position > 44552134-500000 & position < 44580914+500000)
npc1l1 <- vgwas %>% filter(rsid %in% npc1l1$rsid)
npc1l1$ gene <- "NPC1L1"
cetp <- ldl %>% filter(chr == "16" & position > 56995762-500000 & position < 57017757+500000)
cetp <- vgwas %>% filter(rsid %in% cetp$rsid)
cetp$ gene <- "CETP"

# clump mean effects
hmgcr_clump <- ldl %>% filter(rsid %in% hmgcr$rsid) %>% select(rsid, p) %>% rename(pval="p") %>% ld_clump %>% pull(rsid)
hmgcr <- hmgcr %>% filter(rsid %in% hmgcr_clump)
pcsk9_clump <- ldl %>% filter(rsid %in% pcsk9$rsid) %>% select(rsid, p) %>% rename(pval="p") %>% ld_clump %>% pull(rsid)
pcsk9 <- pcsk9 %>% filter(rsid %in% pcsk9_clump)
npc1l1_clump <- ldl %>% filter(rsid %in% npc1l1$rsid) %>% select(rsid, p) %>% rename(pval="p") %>% ld_clump %>% pull(rsid)
npc1l1 <- npc1l1 %>% filter(rsid %in% npc1l1_clump)
cetp_clump <- ldl %>% filter(rsid %in% cetp$rsid) %>% select(rsid, p) %>% rename(pval="p") %>% ld_clump %>% pull(rsid)
cetp <- cetp %>% filter(rsid %in% cetp_clump)

# orient effects to be lipid lowering
ldl <- rbind(hmgcr, pcsk9, npc1l1, cetp)
ldl$flip <- ldl$beta > 0
ldl$allele <- ""
ldl[which(ldl$flip)]$beta <- ldl[which(ldl$flip)]$beta *-1
ldl[which(ldl$flip)]$phi_x <- ldl[which(ldl$flip)]$phi_x *-1
ldl[which(ldl$flip)]$phi_xsq <- ldl[which(ldl$flip)]$phi_xsq *-1
ldl[which(ldl$flip)]$allele <- ldl[which(ldl$flip)]$oa
ldl[which(ldl$flip)]$oa <- ldl[which(ldl$flip)]$ea
ldl[which(ldl$flip)]$ea <- ldl[which(ldl$flip)]$allele
ldl <- ldl %>% select(gene, beta, se, p, phi_x, se_x, phi_xsq, se_xsq, phi_p) %>% mutate_at(vars(beta, se, phi_x, se_x, phi_xsq, se_xsq), funs(round(., 3))) %>% mutate_at(vars(p,phi_p), funs(signif(., 3)))

# HDL-c

# get instruments
hdl <- tophits("ieu-b-109", clump=0)

# load vGWAS for HDL-c
vgwas <- get_variants("hdl_cholesterol.30760.0.0")

# filter drug target loci
cetp <- hdl %>% filter(chr == "16" & position > 56995762-500000 & position < 57017757+500000)
cetp <- vgwas %>% filter(rsid %in% cetp$rsid)
cetp$gene <- "CETP"

# clump mean effects
cetp_clump <- hdl %>% filter(rsid %in% cetp$rsid) %>% select(rsid, p) %>% rename(pval="p") %>% ld_clump %>% pull(rsid)
cetp <- cetp %>% filter(rsid %in% cetp_clump)

# orient effects to be lipid lowering
hdl <- cetp
hdl$flip <- hdl$beta < 0
hdl$allele <- ""
hdl[which(hdl$flip)]$beta <- hdl[which(hdl$flip)]$beta *-1
hdl[which(hdl$flip)]$phi_x <- hdl[which(hdl$flip)]$phi_x *-1
hdl[which(hdl$flip)]$phi_xsq <- hdl[which(hdl$flip)]$phi_xsq *-1
hdl[which(hdl$flip)]$allele <- hdl[which(hdl$flip)]$oa
hdl[which(hdl$flip)]$oa <- hdl[which(hdl$flip)]$ea
hdl[which(hdl$flip)]$ea <- hdl[which(hdl$flip)]$allele
hdl <- hdl %>% select(gene, beta, se, p, phi_x, se_x, phi_xsq, se_xsq, phi_p) %>% mutate_at(vars(beta, se, phi_x, se_x, phi_xsq, se_xsq), funs(round(., 3))) %>% mutate_at(vars(p,phi_p), funs(signif(., 3)))

# write out to file
ldl$trait <- "LDL"
hdl$trait <- "HDL"
all <- rbind(ldl, hdl)
all$phi_1 <- all$phi_x * 1 + all$phi_xsq * (1^2)
all$phi_2 <- all$phi_x * 2 + all$phi_xsq * (2^2)
all$beta.lci <- all$beta - (1.96* all$se)
all$beta.uci <- all$beta + (1.96* all$se)
write.table(all, file="statin.vqtl.txt", sep="\t", quote=F, row.names=F)