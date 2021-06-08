library('data.table')
library('dplyr')
library('purrr')
library('rbgen')

get_filtered_linker <- function(drop_standard_excl=TRUE, drop_non_white_british=TRUE, drop_related=TRUE, application="16729") {
    message("Preparing IEU linker")

    if (application == "16729"){
        # load IEU identifiers to application 16729 linker
        linker=fread("/mnt/storage/private/mrcieu/data/ukbiobank/phenotypic/applications/16729/2019-04-29/data/ieu_id_linker.csv")
        withdrawn=read.table(pipe( 'ssh bc3 "cat /projects/MRC-IEU/research/data/ukbiobank/phenotypic/applications/16729/withdrawals/meta.withdrawn.20200820.csv"' ), stringsAsFactors=F)$V1
    } else if (application == "15825"){
        # load IEU identifiers to application 15825 linker
        linker <- fread("/mnt/storage/private/mrcieu/data/ukbiobank/phenotypic/applications/15825/2019-05-02/data/linker_app15825.csv", col.names=c("appieu", "app15825"))
        withdrawn <- read.table(pipe( 'cat /mnt/storage/private/mrcieu/data/ukbiobank/phenotypic/applications/15825/2019-05-02/data/withdrawals/*' ), stringsAsFactors=F)$V1
    } else {
        abort(paste0("No method for application: ", application))
    }

    # load exclusions
    standard=read.table("/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/derived/standard_exclusions/data.combined_recommended.qctools.txt", stringsAsFactors=F)$V1
    minrelated=read.table("/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/derived/related/relateds_exclusions/data.minimal_relateds.qctools.txt", stringsAsFactors=F)$V1
    highrelated=read.table("/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/derived/related/relateds_exclusions/data.highly_relateds.qctools.txt", stringsAsFactors=F)$V1
    nonwhite=read.table("/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/derived/ancestry/data.non_white_british.qctools.txt", stringsAsFactors=F)$V1

    # exclude samples
    linker=linker[which(!(get(paste0("app", application), linker) %in% withdrawn)),]
    
    if (drop_non_white_british){
        linker=linker[which(!(linker$appieu %in% nonwhite)),]
    }  

    if (drop_standard_excl){
        linker=linker[which(!(linker$appieu %in% standard)),]
    }

    if (drop_related){
        linker=linker[which(!(linker$appieu %in% minrelated)),]
        linker=linker[which(!(linker$appieu %in% highrelated)),]
    }

    message(paste0("Left ", nrow(linker), " samples after exclusions."))

    return(linker)
}

get_genetic_principal_components = function(){
    message("Extracting genetic PCs")
    d = fread("/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/derived/principal_components/data.pca1-40.qctools.txt", header=FALSE, fill = FALSE, sep=" ")
    names(d) = c("appieu", paste0("PC", seq(1, 40)))
    return(d)
}

get_covariates <- function() {
    cov <- fread("/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/derived/standard_covariates/data.covariates.qctools.txt", header=FALSE, fill = FALSE, sep=" ")
    names(cov) = c("appieu", "sex", "chip")
    return(cov)
}

# Extract dosage for variant from BGEN file. Assumes all variants are diploid.
#'
#' @param chrom A single chromosome
#' @param pos Basepair location
#' @param ref Reference allele aka non-effect allele or major allele
#' @param alt Alternative allele aka effect allele or minor allele
#'
#' @return Dataframe of dosages for supplied variant
extract_variant_from_bgen <- function(chrom, pos, ref, alt){
    
    # check inputs
    stopifnot(length(chrom) == 1 && typeof(chrom) == "character")
    stopifnot(length(pos) == 1 && typeof(pos) == "double")
    stopifnot(length(ref) == 1 && typeof(ref) == "character")
    stopifnot(length(alt) == 1 && typeof(alt) == "character")
    # format chrom as expected
    chrom <- gsub("^chr", "", chrom)
    chrom_n <- as.numeric(chrom)
    if (!is.na(chrom_n) && chrom_n < 10){
        chrom <- paste0("0", chrom_n)
    }
    message(paste0("Extracting variant: ", chrom, ":", pos, ref, ">", alt, " from BGEN file"))
    
    # extract variant from bgen
    variant <- bgen.load(paste0("/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr", chrom, ".bgen"),
        data.frame(chromosome = chrom, start = pos, end = pos)
    )
    # check we have a single variant at the correct locus
    if (nrow(variant$variants) == 0){
        stop("Variant not found")
    }
    if (nrow(variant$variants) > 1){
        stop("Variant is multiallelic and not currently supported")
    }
    if (variant$variants$chromosome != chrom){
        stop("Wrong variant returned")
    }
    if (variant$variants$position != pos){
        stop("Wrong variant returned")
    }
    # harmonise
    if (variant$variants$allele0 == ref && variant$variants$allele1 == alt){
        # convert dosage for each genotype to copies of alt allele
        dosage <- as.data.frame(
            apply(variant$data, 1, function(data) { return(data[,1]*0 + data[,2]*1 + data[,3]*2) })
        )
    } else if (variant$variants$allele1 == ref && variant$variants$allele0 == alt){
        # convert dosage for each genotype to copies of alt allele
        dosage <- as.data.frame(
            apply(variant$data, 1, function(data) { return(data[,1]*2 + data[,2]*1 + data[,3]*0) })
        )
    } else {
        stop("Alleles do not match input")
    }
    # tidy up
    names(dosage) <- paste0("chr", gsub("^0", "", chrom), "_", pos, "_", ref, "_", alt)
    dosage$appieu <- row.names(dosage)
    row.names(dosage) <- NULL
    return(dosage)
}

get_variants <- function(trait){
    # load vGWAS & SNP stats; QC loci
    data <- data.frame()

    for (chr in seq(1,22)){
        message(paste0("loading chr", chr))
        if (chr < 10){
            gwas <- fread(paste0("data/", trait, ".vgwas.chr0", chr, ".txt"))
            snp_stats <- fread(paste0("/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/snp-stats/data.chr0", chr, ".snp-stats"), skip=15)
        } else {
            gwas <- fread(paste0("data/", trait, ".vgwas.chr", chr, ".txt"))
            snp_stats <- fread(paste0("/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/snp-stats/data.chr", chr, ".snp-stats"), skip=15)
        }

        # drop multiallelics by rsid
        counts <- table(snp_stats$rsid)
        ma <- as.data.frame(counts[which(counts>1)])
        snp_stats <- snp_stats[!snp_stats$rsid %in% ma$Var1]

        # drop multiallelics by position
        counts <- table(snp_stats$position)
        ma <- as.data.frame(counts[which(counts>1)])
        snp_stats <- snp_stats[!snp_stats$position %in% ma$Var1]
        
        # exclude MAF < 0.05
        snp_stats <- snp_stats[which(snp_stats$minor_allele_frequency > 0.05)]
        
        # exclude HWE violations
        snp_stats <- snp_stats[which(snp_stats$HW_exact_p_value > 1e-5)]

        # exclude high missingness
        snp_stats <- snp_stats[which(snp_stats$missing_proportion < 0.05)]

        # exclude low imputation quality
        snp_stats <- snp_stats[which(snp_stats$info > 0.3)]

        # drop HLA region
        snp_stats <- snp_stats[!(snp_stats$chromosome == 6 & snp_stats$position >= 28477797 & snp_stats$position <= 33448354),]

        # drop vGWAS failed rows
        gwas <- gwas %>% filter(phi_p != -1)

        # merge/filter vGWAS
        snp_stats$key <- paste0(snp_stats$chromosome, "_", snp_stats$position, "_", snp_stats$alleleB, "_", snp_stats$alleleA)
        snp_stats$rsid <- NULL
        gwas$key <- paste0(gwas$chr, "_", gwas$pos, "_", gwas$oa, "_", gwas$ea)
        gwas <- merge(gwas, snp_stats, "key")

        # store results
        data <- rbind(data, gwas)
    }

    return(data)
}

biomarkers <- c(
    "alanine_aminotransferase.30620.0.0",
    "albumin.30600.0.0",
    "alkaline_phosphatase.30610.0.0",
    "apolipoprotein_a.30630.0.0",
    "apolipoprotein_b.30640.0.0",
    "aspartate_aminotransferase.30650.0.0",
    "c_reactive_protein.30710.0.0",
    "calcium.30680.0.0",
    "cholesterol.30690.0.0",
    "creatinine.30700.0.0",
    "cystatin_c.30720.0.0",
    "direct_bilirubin.30660.0.0",
    "gamma_glutamyltransferase.30730.0.0",
    "glucose.30740.0.0",
    "glycated_haemoglobin.30750.0.0",
    "hdl_cholesterol.30760.0.0",
    "igf_1.30770.0.0",
    "ldl_direct.30780.0.0",
    "lipoprotein_a.30790.0.0",
    "oestradiol.30800.0.0",
    "phosphate.30810.0.0",
    "rheumatoid_factor.30820.0.0",
    "shbg.30830.0.0",
    "testosterone.30850.0.0",
    "total_bilirubin.30840.0.0",
    "total_protein.30860.0.0",
    "triglycerides.30870.0.0",
    "urate.30880.0.0",
    "urea.30670.0.0",
    "vitamin_d.30890.0.0"
)

biomarkers_abr <- c(
    "ALT",
    "ALB",
    "ALP",
    "ApoA",
    "ApoB",
    "AST",
    "CRP",
    "Calcium",
    "TC",
    "Creatinine",
    "Cystatin C",
    "Direct bilirubin",
    "GGT",
    "Glucose",
    "HbA1C",
    "HDL",
    "IGF-1",
    "LDL",
    "LipoA",
    "Oestradiol",
    "Phosphate",
    "RF",
    "SHBG",
    "Testosterone",
    "Total bilirubin",
    "Protein",
    "TG",
    "Urate",
    "Urea",
    "Vitamin D"
)
