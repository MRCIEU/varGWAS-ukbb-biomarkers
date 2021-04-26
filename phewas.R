library("data.table")
library("ieugwasr")
library("dplyr")
library('optparse')
library("TwoSampleMR")
library("broom")
set.seed(1234)

option_list = list(
  make_option(c("-t", "--trait"), type="character", default=NULL, help="Name of trait", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# load vGWAS for biomarker risk factor
gwas <- data.frame()
for (chr in seq(1,22)){
    if (chr < 10){
        file <- paste0("data/", opt$trait, ".vgwas.chr0", chr, ".txt")
    } else {
        file <- paste0("data/", opt$trait, ".vgwas.chr", chr, ".txt")
    }
    gwas <- rbind(gwas, fread(file))
}

# load B-P vQTLs
vqtl <- fread(paste0("data/", opt$trait, ".validate.txt"))

# select SNPs which are strongly associated using B-F
vqtl <- vqtl[vqtl$Pvar < 5e-5]

# select vQTLs with evidence of a mean effect on biomaker
mvqtl <- vqtl[vqtl$Pmu < 0.05]

message(paste0("Found ", nrow(mvqtl), " independent mvQTLs"))

# loop over mvQTL
results <- data.frame()
for (snp in mvqtl$rsid){
    message(paste0("Working on ", snp))

    # phewas and take top ten traits
    outcomes <- phewas(snp)

    if (nrow(outcomes) >= 10){
        outcomes <- outcomes[1:10,]
    } else if (nrow(outcomes) == 0){
        next
    }

    # instrument each trait
    for (i in 1:nrow(outcomes)){
        message(paste0("Working on ", outcomes$trait[i]))

        # get instruments for this trait
        iv <- tophits(outcomes$id[i])

        if (nrow(iv) == 0){
            next
        }
        
        # calc absolute Z score
        iv$t <- abs(iv$beta / iv$se)

        # extract vQTLs for instruments
        gwas.outcome <- gwas[gwas$RSID %in% iv$rsid,]

        # add IV-exp effect size
        gwas.outcome <- merge(gwas.outcome, iv[,c("rsid", "t")], by.x="RSID", by.y="rsid")

        # regress IV-outcome F statistic on absolute IV-exp t-score 
        fit <- lm("F ~ t", data=gwas.outcome)

        # add to results
        results <- rbind(results, data.frame(snp, nsnp=nrow(iv), p=tidy(fit)$p.value[2], id=outcomes$id[i], trait=outcomes$trait[i]))
    }
}

print(results[order(results$p, decreasing=F),])