library("data.table")
source("funs.R")
set.seed(23)

# GxE
int <- fread("data/gxe.txt")
int$snp <- stringr::str_split(int$term, ":", simplify=T)[,2]
int$modifier <- stringr::str_split(int$term, ":", simplify=T)[,1]
int$key <- paste0(int$term, ":", int$trait)
int <- int %>% dplyr::filter(p.value < 5e-8)
int_log <- fread("data/gxe-log.txt")
int_log$key <- paste0(int_log$term, ":", int_log$trait)
int_log <- int_log %>% dplyr::filter(p.value < 5e-8)
int <- int %>% dplyr::filter(key %in% int_log$key)

# merge on to gene info
lookup <- fread("Table S1.csv", select=c("snp", "gene", "rsid"))
lookup <- unique(lookup)
int <- merge(int, lookup,"snp")
int$gene <- stringr::str_split(int$gene, "\\|", simplify=T)[,1]
int$modifier <- gsub("age_at_recruitment.21022.0.0", "Age", int$modifier)
int$modifier <- gsub("alcohol_intake_frequency.1558.0.0", "Alcohol", int$modifier)
int$modifier <- gsub("body_mass_index.21001.0.0", "BMI", int$modifier)
int$modifier <- gsub("sex.31.0.0", "Sex", int$modifier)
int$modifier <- gsub("smoking_status.20116.0.0", "Smoking", int$modifier)
int$x <- paste0(int$gene, "(", int$rsid, ")")
int$y <- sapply(int$trait, function(z) biomarkers_abr[z == biomarkers])
int$lci <- int$estimate - (1.96 * int$std.error)
int$uci <- int$estimate + (1.96 * int$std.error)
int <- int %>% dplyr::select(x, modifier, y, estimate, lci, uci, p.value) %>% dplyr::rename(u="modifier")

# GxG
d <- fread("data/gxg.txt")
d$key <- paste0(d$term, ":", d$trait)
d <- d %>% dplyr::filter(p.value < 5e-8)
d_log <- fread("data/gxg-log.txt")
d_log$key <- paste0(d_log$term, ":", d_log$trait)
d_log <- d_log %>% dplyr::filter(p.value < 5e-8)
d <- d %>% dplyr::filter(key %in% d_log$key)

d$lci <- d$estimate - (1.96* d$std.error)
d$uci <- d$estimate + (1.96* d$std.error)
d <- d %>% dplyr::filter(p.value < 5e-8)
d <- cbind(d, data.frame(stringr::str_split(d$term, ":", simplify=T), stringsAsFactors=F))
d <- merge(d, lookup, by.x="X1", by.y="snp")
d <- merge(d, lookup, by.x="X2", by.y="snp")
d$x <- paste0(d$rsid.x, "(", d$gene.x, ")")
d$u <- paste0(d$rsid.y, "(", d$gene.y, ")")
d$y <- sapply(d$trait, function(z) biomarkers_abr[z == biomarkers])
d <- d %>% dplyr::select(x, u, y, estimate, lci, uci, p.value) %>% dplyr::rename()

# combine
all <- rbind(int, d)
names(all) <- c(
    "snp",
    "modifier",
    "outcome",
    "beta",
    "lci",
    "uci",
    "p"
)

# save
write.csv(all, file="Table S2.csv", row.names=F, quote=F)