library("data.table")
library("ggplot2")
library("dplyr")
library("stringr")
library("ggpubr")
library("viridis")
library("lmtest")
library("sandwich")
library('forestplot')
library("viridis")
library("RColorBrewer")
library("grid")
library("broom")
library("gwasglue")
library("ieugwasr")
library("susieR")
library("robustbase")
source("funs.R")
options(ieugwasr_api="http://64.227.44.193:8006/")
set.seed(123)

get_dat <- function(file){
    # read in gxe results
    d <- fread(file)

    # drop BMI
    d <- d %>% dplyr::filter(trait != "body_mass_index.21001.0.0")

    # filter SNPs to show
    d <- d %>% dplyr::filter(p.value < 5e-8)

    # merge
    d <- cbind(d, as.data.frame(stringr::str_split(d$term, ":", simplify=T), stringsAsFactors=F))
    
    # add key
    d$tt <- paste0(d$trait, ":", d$term)

    return(d)
}

finemap_func <- function(chr_pos, id){
    message(paste0("Working on variant: ", chr_pos))

    # select natural interval around lead SNP
    region <- map_variants_to_regions(chrpos=chr_pos, pop="EUR")

    # select variants and LD matrix for region
    dat <- ieugwasr_to_finemapr(
        region$region,
        id,
        bfile = "/mnt/storage/home/ml18692/projects/varGWAS-ukbb-biomarkers/data/EUR",
        plink_bin = "/mnt/storage/home/ml18692/projects/varGWAS-ukbb-biomarkers/data/plink_Linux"
    )

    # perform finemapping using SuSie
    fitted_rss <- tryCatch(
        expr = {
            susieR::susie_rss(
                dat[[1]]$z$zscore,
                dat[[1]]$ld,
                L=10,
                estimate_prior_variance=TRUE
            )
        },
        error = function(e){ 
            NULL
        }
    )

    if (is.null(fitted_rss)){
        return(NULL)
    }

    # collect fine mapped snps
    cs <- summary(fitted_rss)$cs

    if (is.null(cs)){
        return(NULL)
    }

    snps <- data.frame()
    for (j in 1:nrow(cs)){
        for (k in stringr::str_split(cs$variable[j], ",", simplify=T) %>% as.numeric){
            res <- cs[j,]
            res <- cbind(res, data.frame(
                snp=dat[[1]]$z$snp[k],
                pip=fitted_rss$pip[k]
            ))
            snps <- rbind(snps, res)
        }
    }

    # get assoc for snps
    ass <- associations(snps$snp, id) %>% dplyr::select(rsid, chr, position, nea, ea) %>% dplyr::rename(oa="nea", snp="rsid")
    snps <- merge(snps, ass, "snp")

    return(data.frame(chr_pos, snps, id, stringsAsFactors=F))
}

get_finemap <- function(d){
    # add opengwas ID
    d$id <- paste0("ukb-d-", str_split(d$trait, "\\.", simplify=T)[,2], "_irnt")

    # extract chr-pos
    V1 <- as.data.frame(str_split(d$V1, "_", simplify=T), stringsAsFactors=F)
    names(V1) <- c("chr.1", "pos.1", "oa.1", "ea.1")
    V1$chr.1 <- gsub("chr", "", V1$chr.1)
    V2 <- as.data.frame(str_split(d$V2, "_", simplify=T), stringsAsFactors=F)
    names(V2) <- c("chr.2", "pos.2", "oa.2", "ea.2")
    V2$chr.2 <- gsub("chr", "", V2$chr.2)    
    d <- cbind(d, V1, V2)
    d$chr_pos.1 <- paste0(d$chr.1, ":", d$pos.1)
    d$chr_pos.2 <- paste0(d$chr.2, ":", d$pos.2)

    results <- data.frame()
    for (i in 1:nrow(d)){
        res <- finemap_func(d$chr_pos.1[i], d$id[i])
        if (!is.null(res)){
            res$trait <- d$trait[i]
            results <- rbind(res, results)
        }
        res <- finemap_func(d$chr_pos.2[i], d$id[i])
        if (!is.null(res)){
            res$trait <- d$trait[i]
            results <- rbind(res, results)
        }
    }

    return(results)
}

get_est <- function(trait, v1, v2, finemap){
    # read in extracted phenotypes
    pheno <- fread(paste0("data/", trait, ".txt"))

    # subset finemapped variants for this trait
    finemap <- finemap %>% dplyr::filter(trait == !!trait)

    # load dosages
    snps1 <- stringr::str_split(v1, "_", simplify=T) %>% as.data.frame(., stringsAsFactors=F) %>% dplyr::rename(chr="V1", pos="V2", oa="V3", ea="V4") %>% dplyr::mutate(chr=gsub("chr", "", chr))
    snps2 <- stringr::str_split(v2, "_", simplify=T) %>% as.data.frame(., stringsAsFactors=F) %>% dplyr::rename(chr="V1", pos="V2", oa="V3", ea="V4") %>% dplyr::mutate(chr=gsub("chr", "", chr))
    snps3 <- finemap %>% dplyr::select(chr, position, oa, ea) %>% dplyr::rename(pos="position")
    snps <- unique(rbind(snps1, snps2, snps3))
    snps$key <- paste0("chr", snps$chr, "_", snps$pos, "_", snps$oa, "_", snps$ea)
    chr_pos.1 <- paste0(snps1$chr, ":", snps1$pos)
    chr_pos.2 <- paste0(snps2$chr, ":", snps2$pos)

    for (i in 1:nrow(snps)){
        dosage <- extract_variant_from_bgen(as.character(snps$chr[i]), as.double(snps$pos[i]), snps$oa[i], snps$ea[i])
        pheno <- merge(pheno, dosage, "appieu")
    }

    # test for interaction adjusting for finemapped variants
    adj1 <- finemap %>% dplyr::filter(chr_pos == !!chr_pos.1) %>% dplyr::select(chr, position, oa, ea) %>% dplyr::mutate(snp=paste0("chr", chr, "_", position, "_", oa, "_", ea)) %>% pull(snp)
    adj2 <- finemap %>% dplyr::filter(chr_pos == !!chr_pos.2) %>% dplyr::select(chr, position, oa, ea) %>% dplyr::mutate(snp=paste0("chr", chr, "_", position, "_", oa, "_", ea)) %>% pull(snp)
    adj <- unique(adj1, adj2)
    f <- paste0(trait, " ~ ", v1, " * ", v2, " + age_at_recruitment.21022.0.0 + sex.31.0.0 + ", paste0("PC", seq(1,10), collapse=" + "), " + ", paste0(adj, collapse=" + "))
    message(f)
    mod <- lm(f, data=pheno)
    t <- coeftest(mod, vcov = vcovHC(mod, type = "HC0")) %>% tidy
    t <- t %>% dplyr::filter(grepl(":", t$term))
    t$trait <- trait
    t$formula <- f
    return(t)
}

get_plot <- function(d){
    # add rsid
    lookup <- data.frame(
        term=c(
            "chr22_44324730_T_C:chr4_88212722_A_G",
            "chr19_45411941_C_T:chr12_121420260_G_A",
            "chr19_19379549_T_C:chr19_45416741_T_C",
            "chr9_136149830_A_G:chr19_49210869_G_A",
            "chr11_116648917_C_G:chr8_19863507_T_C",
            "chr8_19912370_A_G:chr11_116651115_T_C"
        ), 
        f=c(
            "PNPLA3 (rs738408C) x HSD17B13 (rs13141441G)",
            "APOE (rs429358T) x HNF1A (rs7979473A)",
            "TM6SF2 (rs58542926C) x APOE (rs438811C)",
            "ABO (rs532436G) x FUT2 (rs2638281A)",
            "APOA5 (rs964184) x LPL (rs2119689C)",
            "LPL (rs115849089G) x APOA5 (rs11604424C)"
        )
    )
    d <- merge(d, lookup, "term", all.T=T)

    d$Trait <- sapply(d$trait, function(x) return(biomarkers_abr[biomarkers==x]))

    # create row key
    key <- data.frame(Trait=sort(unique(d$Trait)), stringsAsFactors=F)
    key$key <- row(key) %% 2
    d <- merge(d, key, "Trait")
    d$key <- factor(d$key)

    # count number of traits
    n_traits <- length(unique(d$Trait))

    # Create a data frame with the faceting variables
    # and some dummy data (that will be overwritten)
    tp <- data.frame()
    for (tr in unique(d$Trait)){
        tp <- rbind(tp, data.frame(
            Trait=tr,
            fill=which(tr == unique(d$Trait)) %% 2
        ))
    }
    tp$fill <- as.factor(tp$fill)

    # create plot
    p <- ggplot(d, aes(x=f, y=est, ymin=lci, ymax=uci)) +
        coord_flip() +
        facet_grid(Trait~., scales="free", space="free_y") +
        geom_point(size = 1.5) +
        geom_errorbar(width=.05) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
        theme_classic() +
        scale_y_continuous(limits = c(-.1, .1), breaks=c(-.1, -0.05, 0, 0.05, .1)) +
        geom_rect(inherit.aes = F, show.legend = FALSE, data = tp, aes(fill = fill), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.15) +
        scale_fill_manual(values=brewer.pal(2,"Paired")) +
        theme(
            axis.title.y = element_blank(),
            strip.background = element_blank(),
            strip.text.y = element_text(angle = 0),
            legend.position = "bottom",
            panel.spacing.y = unit(0, "lines")
        ) +
        ylab("Genotype * genotype (dosage) interaction effect estimate, SD (95% CI)")
    return(p)
}

# load gxe effects
additive <- get_dat("data/gxg.txt")
multiplicative <- get_dat("data/gxg-log.txt")

# append sensitivity P value
additive <- merge(additive, fread("data/gxg-log.txt") %>% dplyr::mutate(tt=paste0(trait, ":", term)) %>% dplyr::select(tt, p.value) %>% dplyr::rename(p_sens="p.value"), "tt")
multiplicative <- merge(multiplicative, fread("data/gxg.txt") %>% dplyr::mutate(tt=paste0(trait, ":", term)) %>% dplyr::select(tt, p.value) %>% dplyr::rename(p_sens="p.value"), "tt")

# finemap gxg snps
additive.finemap <- get_finemap(additive)
multiplicative.finemap <- get_finemap(multiplicative)

# test for effect adjusting for finemapped variants
additive.results <- data.frame()
for (i in 1:nrow(additive)){
    est <- get_est(additive$trait[i], additive$V1[i], additive$V2[i], additive.finemap)
    additive.results <- rbind(additive.results, est)
}
multiplicative.results <- data.frame()
for (i in 1:nrow(multiplicative)){
    est <- get_est(multiplicative$trait[i], multiplicative$V1[i], multiplicative$V2[i], multiplicative.finemap, multiplicative=T)
    multiplicative.results <- rbind(multiplicative.results, est)
}

# plot
pdf("gxg-additive-finemap.pdf", height=6, width=8)
print(get_plot(additive.results))
dev.off()

write.table(additive.results, sep="\t", quote=F, row.names=F, file=paste0("data/", opt$trait, ".gxg-add-finemap.txt"))

pdf("gxg-multiplicative-finemap.pdf", height=5, width=8)
print(get_plot(multiplicative.results))
dev.off()

write.table(multiplicative.results, sep="\t", quote=F, row.names=F, file=paste0("data/", opt$trait, ".gxg-multi-finemap.txt"))