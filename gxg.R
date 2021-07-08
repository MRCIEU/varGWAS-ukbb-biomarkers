library('optparse')
library('data.table')
library('dplyr')
library('broom')
library('ieugwasr')
library('robustbase')
source("funs.R")
set.seed(1234)
options(ieugwasr_api="http://64.227.44.193:8006/")

option_list = list(
  make_option(c("-t", "--trait"), type="character", default=NULL, help="Name of trait", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

message(paste0("trait ", opt$trait))

# read in extracted phenotypes
pheno <- fread(paste0("data/", opt$trait, ".txt"))

# read in clumped vQTLs
snps <- fread(paste0("data/", opt$trait, ".clump.txt"))

# add key
snps$key <- paste0("chr", snps$chr, "_", snps$pos, "_", snps$oa, "_", snps$ea)

# load dosages
for (i in 1:nrow(snps)){
    dosage <- extract_variant_from_bgen(as.character(snps$chr[i]), as.double(snps$pos[i]), snps$oa[i], snps$ea[i])
    pheno <- merge(pheno, dosage, "appieu")
}

# select vQTLs
vqtls <- grep("^chr", names(pheno), value=T)

# test for interaction between each snp
results <- data.frame()
for (i in 1:length(vqtls)){
  for (j in 1:length(vqtls)){

    # skip GxG on same chromosome within 10Mb
    i_chr <- snps %>% filter(key == vqtls[i]) %>% pull("chr")
    i_pos <- snps %>% filter(key == vqtls[i]) %>% pull("pos")
    j_chr <- snps %>% filter(key == vqtls[j]) %>% pull("chr")
    j_pos <- snps %>% filter(key == vqtls[j]) %>% pull("pos")

    if (i_chr == j_chr){
      if (abs(i_pos - j_pos) < 10000000){
        message("Skipping test (<10Mb) for: ", vqtls[i], " ", vqtls[j])
        next
      }
    }

    if (paste0(vqtls[j], ":" ,vqtls[i]) %in% results$term){
        message("Skipping test (already done) for: ", vqtls[i], " ", vqtls[j])
        next
    }

    # test GxG
    message("Testing GxG for: ", vqtls[i], " ", vqtls[j])
    f <- as.formula(paste0(opt$trait, " ~ age_at_recruitment.21022.0.0 + sex.31.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + ", paste0(vqtls[i], " * " ,vqtls[j], collapse=" + ")))
    t <- tidy(lmrob(f, data=pheno))

    # store results
    results <- rbind(results, t[grep(":", t$term),])
  }
}

# save
write.table(results, sep="\t", quote=F, row.names=F, file=paste0("data/", opt$trait, ".gxg.txt"))