library("data.table")

# load in SNP estimates
r <- fread("data/ukb_bmi.vgwas.r_subsample.txt")
cpp <- fread("data/ukb_bmi.vgwas.chr22.txt")

# merge data
cpp$key <- paste0("chr", cpp$CHR, "_", cpp$POS, "_", cpp$OA, "_", cpp$EA)
d <- merge(r, cpp, by.x="term", by.y="key")

# check correlation between R and CPP
cor.test(d$BETA, d$estimate)
cor.test(d$SE, d$std.error)
cor.test(d$P, d$p.value)