library("data.table")
library("dplyr")
set.seed(1234)

results <- data.frame()

# use exp IVs and estimate variance effect with combination for 1-or-2 copies of allele
# forest plot mean and variance effects
# dicuss issues around gamma distributions (mean-variance rel) & variance redcuing effect which is good vs increasing variance that wud suggest subgroups & reference Senn & cortes papers

# load in clumped vQTLs
d <- fread("data/vqtls.txt")

# select within .5Mb of lipid lowering targets
res <- d %>% dplyr::filter(chr == 5 & position > 74632154-500000 & position < 74657929+500000)
res$target <- "HMGCR"
results <- rbind(results, res)

res <- d %>% dplyr::filter(chr == 1 & position > 55505221-500000 & position < 55530525+500000)
res$target <- "PCSK9"
results <- rbind(results, res)

res <- d %>% dplyr::filter(chr == 7 & position > 44552134-500000 & position < 44580914+500000)
res$target <- "NPC1L1"
results <- rbind(results, res)

res <- d %>% dplyr::filter(chr == 16 & position > 56995762-500000 & position < 57017757+500000)
res$target <- "CETP"
results <- rbind(results, res)

results %>% select(trait, rsid, target, phi_p)