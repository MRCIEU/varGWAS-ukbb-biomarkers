library('optparse')
library('data.table')
library('dplyr')
library('broom')
library('ieugwasr')
library("lmtest")
library("sandwich")
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

# read in main clumped vQTLs P < 5e-8/30 & add key
main <- fread(paste0("data/", opt$trait, ".clump.txt"))
main$key <- paste0("chr", main$chr, "_", main$pos, "_", main$oa, "_", main$ea)
stopifnot(nrow(main) > 0)

# read in modifier clumped vQTLs P < 5e-5 & add key
modifier <- fread(paste0("data/", opt$trait, ".clump.txt"))
modifier$key <- paste0("chr", modifier$chr, "_", modifier$pos, "_", modifier$oa, "_", modifier$ea)

# load dosages
snps <- main %>% dplyr::select(chr, pos, oa, ea)
snps <- rbind(snps, modifier %>% dplyr::select(chr, pos, oa, ea))
snps <- unique(snps)
snps$key <- paste0("chr", snps$chr, "_", snps$pos, "_", snps$oa, "_", snps$ea)

for (i in 1:nrow(snps)){
    dosage <- extract_variant_from_bgen(as.character(snps$chr[i]), as.double(snps$pos[i]), snps$oa[i], snps$ea[i])
    pheno <- merge(pheno, dosage, "appieu")
}
pheno[[paste0(opt$trait, "_log")]] <- log(pheno[[opt$trait]])

# test for interaction between each snp
results <- data.frame()
results_log <- data.frame()
for (i in 1:length(main$key)){
  for (j in 1:length(modifier$key)){

    # skip GxG on same chromosome within 10Mb
    i_chr <- snps %>% dplyr::filter(key == main$key[i]) %>% dplyr::pull("chr")
    i_pos <- snps %>% dplyr::filter(key == main$key[i]) %>% dplyr::pull("pos")
    j_chr <- snps %>% dplyr::filter(key == modifier$key[j]) %>% dplyr::pull("chr")
    j_pos <- snps %>% dplyr::filter(key == modifier$key[j]) %>% dplyr::pull("pos")

    if (i_chr == j_chr){
      if (abs(i_pos - j_pos) < 10000000){
        message("Skipping test (<10Mb) for: ", main$key[i], " ", modifier$key[j])
        next
      }
    }

    if (paste0(modifier$key[j], ":" ,main$key[i]) %in% results$term){
        message("Skipping test (already done) for: ", main$key[i], " ", modifier$key[j])
        next
    }

    message("Testing GxG for: ", main$key[i], " ", modifier$key[j])

    # test GxG on additive scale
    f <- as.formula(paste0(opt$trait, " ~ age_at_recruitment.21022.0.0 + sex.31.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + ", paste0(main$key[i], " * " ,modifier$key[j], collapse=" + ")))
    mod <- lm(f, data=pheno)
    t <- coeftest(mod, vcov = vcovHC(mod, type = "HC0")) %>% tidy
    t <- t %>% dplyr::filter(grepl(":", t$term))
    results <- rbind(results, t)

    # test GxG on multiplicative scale
    f <- as.formula(paste0(opt$trait, "_log ~ age_at_recruitment.21022.0.0 + sex.31.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + ", paste0(main$key[i], " * " ,modifier$key[j], collapse=" + ")))
    mod <- lm(f, data=pheno)
    t <- coeftest(mod, vcov = vcovHC(mod, type = "HC0")) %>% tidy
    t <- t %>% dplyr::filter(grepl(":", t$term))
    results_log <- rbind(results_log, t)
    
  }
}

# save
write.table(results, sep="\t", quote=F, row.names=F, file=paste0("data/", opt$trait, ".gxg.txt"))
write.table(results_log, sep="\t", quote=F, row.names=F, file=paste0("data/", opt$trait, ".gxg-log.txt"))