library('optparse')
library('data.table')
library('dplyr')
library('broom')
library('ieugwasr')
source("funs.R")
set.seed(1234)
options(ieugwasr_api="http://64.227.44.193:8006/")

get_gxg <- function(id, chra, starta, enda, chrb, startb, endb, trait){
  # get tophits
  hits <- tophits(id=id, pval = 5e-05)
  gene1 <- hits %>% filter(chr == chra & position >= starta-500000 & position <= enda+500000)
  gene2 <- hits %>% filter(chr == chrb & position >= startb-500000 & position <= endb+500000)
  print(gene2$rsid)

  # read in extracted phenotypes
  opt <- data.frame(trait=trait, stringsAsFactors=F)
  pheno <- fread(paste0("data/", opt$trait, ".txt"))

  # add key
  snps <- rbind(gene1, gene2)
  snps$pos <- snps$position
  snps$key <- paste0("chr", snps$chr, "_", snps$pos, "_", snps$nea, "_", snps$ea)

  # load dosages
  for (i in 1:nrow(snps)){
      dosage <- extract_variant_from_bgen(as.character(snps$chr[i]), as.double(snps$pos[i]), snps$nea[i], snps$ea[i])
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
      j_chr <- snps %>% filter(key == vqtls[j]) %>% pull("chr")

      if (i_chr == j_chr){
        next
      }

      if (paste0(vqtls[j], ":" ,vqtls[i]) %in% results$term){
          message("Skipping test (already done) for: ", vqtls[i], " ", vqtls[j])
          next
      }

      # test GxG
      message("Testing GxG for: ", vqtls[i], " ", vqtls[j])
      f <- as.formula(paste0(opt$trait, " ~ age_at_recruitment.21022.0.0 + sex.31.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + ", paste0(vqtls[i], " * " ,vqtls[j], collapse=" + ")))
      fit <- lm(f, pheno)
      t <- tidy(fit)

      # store results
      results <- rbind(results, t[grep(":", t$term),])
    }
  }
  return(results)
}

# test for gxg
tg <- get_gxg("ukb-d-30870_irnt", 11, 116660083, 116663136, 8, 19759228, 19824769, "triglycerides.30870.0.0")
crp <- get_gxg("ukb-d-30710_irnt", 12, 121416346, 121440315, 19, 45409011, 45412650, "c_reactive_protein.30710.0.0")
alt <- get_gxg("ukb-d-30620_irnt", 22, 44319619, 44360368, 4, 88224941, 88244058, "alanine_aminotransferase.30620.0.0")