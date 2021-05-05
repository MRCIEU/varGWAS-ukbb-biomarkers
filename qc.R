library('data.table')
library('dplyr')
library('qqman')
library('stringr')
library('optparse')
set.seed(1234)

option_list = list(
  make_option(c("-t", "--trait"), type="character", default=NULL, help="Name of trait", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

trait_name <- opt$trait
trait_name <- str_split(trait_name, "\\.", simplify = TRUE)[,1]
trait_name <- gsub("_", " ", trait_name)
trait_name <- str_to_title(trait_name)
message(paste0("trait ", opt$trait))
message(paste0("trait name ", trait_name))

# load vGWAS & SNP stats; QC loci
data <- data.frame()
for (chr in seq(1,22)){
    message(paste0("loading chr", chr))
    if (chr < 10){
        gwas <- fread(paste0("data/", opt$trait, ".vgwas.chr0", chr, ".txt"))
        snp_stats <- fread(paste0("/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/snp-stats/data.chr0", chr, ".snp-stats"), skip=15)
    } else {
        gwas <- fread(paste0("data/", opt$trait, ".vgwas.chr", chr, ".txt"))
        snp_stats <- fread(paste0("/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/snp-stats/data.chr", chr, ".snp-stats"), skip=15)
    }
    
    # exclude MAF < 0.05
    snp_stats <- snp_stats[which(snp_stats$minor_allele_frequency > 0.05)]
    
    # exclude HWE violations
    snp_stats <- snp_stats[which(snp_stats$HW_exact_p_value > 1e-5)]

    # exclude high missingness
    snp_stats <- snp_stats[which(snp_stats$missing_proportion < 0.05)]

    # exclude low imputation quality
    snp_stats <- snp_stats[which(snp_stats$info > 0.3)]

    # drop multiallelics by rsid
    counts <- table(snp_stats$rsid)
    ma <- as.data.frame(counts[which(counts>1)])
    snp_stats <- snp_stats[!snp_stats$rsid %in% ma$Var1]

    # drop multiallelics by position
    counts <- table(snp_stats$position)
    ma <- as.data.frame(counts[which(counts>1)])
    snp_stats <- snp_stats[!snp_stats$position %in% ma$Var1]

    # drop HLA region
    snp_stats <- snp_stats[!(snp_stats$chromosome == 6 & snp_stats$position >= 28477797 & snp_stats$position <= 33448354),]

    # drop vGWAS failed rows
    gwas <- gwas %>% filter(P != -1)

    # merge/filter vGWAS
    snp_stats$key <- paste0(snp_stats$chromosome, "_", snp_stats$position, "_", snp_stats$alleleA, "_", snp_stats$alleleB)
    gwas$key <- paste0(gwas$CHR, "_", gwas$POS, "_", gwas$OA, "_", gwas$EA)
    gwas <- merge(gwas, snp_stats, "key")

    # store results
    data <- rbind(data, gwas)
}

# manhattan
png(paste0("data/", opt$trait, "_phi_manhattan.png"))
manhattan(data, ylim = c(0, 25), chr="CHR", bp="POS", p="P", snp="RSID", main = trait_name)
dev.off()

# qq plot
png(paste0("data/", opt$trait, "_phi_qq.png"))
qq(data$P, main = trait_name)
dev.off()