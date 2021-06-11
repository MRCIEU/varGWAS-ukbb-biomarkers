load("data/pheno.RData")
library('optparse')
library('data.table')
library('dplyr')
library('broom')
library("stringr")
library('forestplot')
source("funs.R")
set.seed(1234)

option_list = list(
  make_option(c("-t", "--trait"), type="character", default=NULL, help="Name of trait", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

message(paste0("trait ", opt$trait))

# load linker
linker <- get_filtered_linker(drop_standard_excl=TRUE, drop_non_white_british=FALSE, drop_related=FALSE, application="15825")

# merge data
dat <- merge(linker, pheno, by.x="app15825", by.y="eid")

# keep sens cols
dat <- dat %>% select(appieu, place_of_birth_in_UK_north_co_ordinate.129.0.0, place_of_birth_in_UK_east_co_ordinate.130.0.0)

# read in extracted phenotypes
pheno <- fread(paste0("data/", opt$trait, ".txt"))
pheno <- merge(pheno, dat, "appieu")

# read in main GxG associations
snps <- fread(paste0("data/", opt$trait, ".gxg.txt"))

# take top n=20 GxG hits for replication
snps <- snps[order(snps$p.value)] %>% head(n=20)

# split term
snps <- cbind(snps, str_split(snps$term, ":", simplify=T), stringsAsFactors=F)
usnps <- unique(c(snps$V1,snps$V2), stringsAsFactors=F)
usnps <- as.data.frame(str_split(usnps, "_", simplify=T), stringsAsFactors=F)
usnps$V1 <- gsub("chr", "", usnps$V1)
usnps$V2 <- as.numeric(usnps$V2)

# load dosages
for (i in 1:nrow(usnps)){
    dosage <- extract_variant_from_bgen(usnps$V1[i], usnps$V2[i], usnps$V3[i], usnps$V4[i])
    pheno <- merge(pheno, dosage, "appieu")
}

# test for interaction between each snp
results <- data.frame()
for (i in 1:length(snps$term)){
  # test GxG
  pair <- str_split(snps$term[i], ":", simplify=T)
  message("Testing GxG for: ", pair[1], " ", pair[2])
  f <- as.formula(paste0(opt$trait, " ~ age_at_recruitment.21022.0.0 + sex.31.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + chip + place_of_birth_in_UK_north_co_ordinate.129.0.0 + place_of_birth_in_UK_east_co_ordinate.130.0.0 + ", paste0(pair[1], " * " ,pair[2], collapse=" + ")))
  fit <- lm(f, pheno)
  t <- tidy(fit)
  
  # store results
  results <- rbind(results, t[grep(":", t$term),])
}

# load within fam analysis
wf <- fread(paste0("data/", opt$trait, ".gxg-wf.txt"))

# combine with main analysis
results <- merge(snps, results, "term", suffixes=c(".main", ".sens"))
results <- merge(results, wf, "term")
results$estimate.main.lci <- results$estimate.main - (1.96*results$std.error.main)
results$estimate.main.uci <- results$estimate.main + (1.96*results$std.error.main)
results$estimate.sens.lci <- results$estimate.sens - (1.96*results$std.error.sens)
results$estimate.sens.uci <- results$estimate.sens + (1.96*results$std.error.sens)
results$beta.lci <- results$beta - (1.96*results$se)
results$beta.uci <- results$beta + (1.96*results$se)

# forest plot
forestplot(
    results$term,
    boxsize = 0.05,
    xticks = c(-.5, -.25, 0, .25, .5),
    legend = c("Main", "Senstivity", "Within-family"),
    mean = results[,c("estimate.main", "estimate.sens", "beta")],
    lower = results[,c("estimate.main.lci","estimate.sens.lci", "beta.lci")],
    upper = results[,c("estimate.main.uci", "estimate.sens.uci", "beta.uci")],
    col=fpColors(lines=c("#3182bd", "#a28bf7", "#31a354"), box=c("#3182bd","#a28bf7", "#31a354")),
    xlab=paste0("GxG effect on ",opt$trait," (SD [95% CI])"),
    txt_gp = fpTxtGp(ticks=gpar(cex=1), xlab=gpar(cex=1))
)