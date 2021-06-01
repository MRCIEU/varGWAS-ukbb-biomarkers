library("data.table")
library("dplyr")
library("stringr")
library("GenomicRanges")
library("IRanges")
set.seed(123)

# read in MR results
d <- fread("data/mr.txt")

# add key
d$key <- paste0(d$exposure, ".", d$outcome)

# Steiger filter
d <- d[d$correct_causal_direction]

# drop null MR estimates for IVW / Wald
main <- d %>% filter(method == "Inverse variance weighted" | method == "Wald ratio")
main <- main %>% filter(pval < (0.05/nrow(main)))

# drop estimates not supported by weighted median
sens <- d %>% filter(method == "Weighted median")
sens <- sens %>% filter(pval < 0.05)
sig <- main[main$key %in% sens$key | main$method == "Wald ratio"]

# filter out bidirectional effects
bi <- sig[grepl("irnt", sig$exposure)]
sig <- sig[!grepl("irnt", sig$exposure)]
bi$id <- gsub("_irnt", "", gsub("ukb-d-", "", bi$id.exposure))
for (i in 1:nrow(bi)){
    sig <- sig[!(sig$id.exposure == bi$id.outcome[i] & grepl(bi$id[i], sig$outcome))]
}

# add target name
sig$target <- str_split(sig$exposure, " \\|\\| ", simplify=T)[,1]

# load mapping to ENSG
lookup <- read.table(file="data/lookup.grch37.txt", stringsAsFactors=F)

# merge data
sig <- merge(sig, lookup, by.x="target", by.y="Target")
sig$startpad <- sig$start - 500000
sig$endpad <- sig$end + 500000
sig$region <- paste0(sig$chrom, ":", sig$startpad, "-", sig$endpad)

# read in coloc data
coloc <- fread("data/coloc.txt")
coloc <- coloc[coloc$PP.H4.abf > .8]

# loop over MR estimates for effect of gene on biomarker conc
results <- data.frame()
for (i in 1:nrow(sig)){
    # grab coloc vQTLs for cis region of this gene
    ctemp <- coloc %>% filter(trait==!!sig$outcome[i])
    cregions <- GRanges(ctemp$region)
    
    # test for intersection between MR and coloc
    mregion <- GRanges(sig$region[i])

    # count overlaps between coloc and MR data
    counts <- suppressWarnings(countOverlaps(cregions, mregion))

    # save res
    if (sum(counts) > 0){
        res <- ctemp[counts>0]
        res <- cbind(res, sig[i], stringsAsFactors=F)
        results <- rbind(results, res)
    }
}