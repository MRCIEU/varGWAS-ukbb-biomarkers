source("funs.R")
set.seed(1234)

# load clumped vQTLs
results <- fread("data/vqtls.txt")
results$key <- paste0(results$chr, ":", results$position)
results$snp <- paste0("chr", results$chr, "_", results$pos, "_", results$oa, "_", results$ea)

# get data on nearest gene
ng <- fread("data/nearest.txt")
ng <- unique(ng)
ng$key <- paste0(ng$V1,":",ng$V2)
ng$V1 <- NULL
ng$V2 <- NULL
names(ng)[1] <- "gene"
ng <- ng %>%
  group_by_at(vars(key)) %>%
  summarize(gene = toString(gene)) %>%
  ungroup()
ng$gene <- gsub(", ", "|", ng$gene)

# results with nearest gene
all <- merge(results, ng, "key", all.x=T)
all$key <- NULL

# tidy outcome name
all$outcome <- sapply(all$trait, function(x) biomarkers_abr[x==biomarkers], simplify=T)

# select fields for table
all <- all %>% dplyr::select(snp,rsid,outcome,gene,phi_x1,se_x1,phi_x2,se_x2,phi_f,phi_p)

# write to table
write.csv(all, file="Table S1.csv", quote=F, row.names=F)