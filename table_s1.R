source("funs.R")
set.seed(1234)

# load clumped vQTLs
results <- fread("data/vqtls.txt")

# get data on nearest gene
ng <- fread("data/nearest.txt", col.names=c("rsid", "gene"))

# results with nearest gene
all <- merge(results, ng, "rsid", all.x=T)
all$key <- NULL

# tidy outcome name
all$outcome <- sapply(all$trait, function(x) biomarkers_abr[x==biomarkers], simplify=T)

# select fields for table
all <- all %>% dplyr::select(rsid,outcome,gene,phi_x1,se_x1,phi_x2,se_x2,phi_f,phi_p)

# write to table
write.csv(all, file="Table S1.csv", quote=F, row.names=F)