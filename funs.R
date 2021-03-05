library('data.table')
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