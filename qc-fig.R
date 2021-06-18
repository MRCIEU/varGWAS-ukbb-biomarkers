library('data.table')
library('dplyr')
library('stringr')
library('ggpubr')
source("funs.R")
set.seed(1234)

# Taken from https://danielroelfs.com/blog/how-i-create-manhattan-plots-using-ggplot/
man_plot <- function(gwas_data, sig=5e-8/30, ylim=30){
    data_cum <- gwas_data %>% 
        group_by(chr) %>% 
        summarise(max_bp = max(pos)) %>% 
        mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
        select(chr, bp_add)

    gwas_data <- gwas_data %>% 
        inner_join(data_cum, by = "chr") %>% 
        mutate(bp_cum = pos + bp_add)

    axis_set <- gwas_data %>% 
        group_by(chr) %>% 
        summarize(center = mean(bp_cum))

    gwas_data$p <- -log10(gwas_data$p)

    manhplot <- ggplot(gwas_data, aes(x = bp_cum, y = p, color = as.factor(chr), size = p)) +
        geom_point(alpha = 0.75) +
        geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") + 
        scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
        scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
        scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$chr)))) +
        scale_size_continuous(range = c(0.5,3)) +
        labs(x = NULL, y = expression(paste("-log"[10],"(", plain(P),")"))) + 
        theme_minimal() +
        theme( 
            legend.position = "none",
            strip.text.x = element_text(size = 40),
            strip.text = element_text(size = 40),
            panel.border = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            axis.text = element_text(size = 40),
            axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.y = element_text(size = 40),
            axis.title = element_text(size = 40)
        ) +
        facet_wrap(~trait)

    return(manhplot)
}

# Taken from https://danielroelfs.com/blog/how-i-make-qq-plots-using-ggplot/
qq_plot_dat <- function(sumstats.data, pval_col, outcome, ci=0.95){
    nSNPs <- nrow(sumstats.data)
    plotdata <- data.frame(
        observed = -log10(sort(sumstats.data[[pval_col]])),
        expected = -log10(ppoints(nSNPs)),
        clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = seq(nSNPs), shape2 = rev(seq(nSNPs)))),
        cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = seq(nSNPs), shape2 = rev(seq(nSNPs)))),
        outcome
    )
    return(plotdata)
}

# Taken from https://danielroelfs.com/blog/how-i-make-qq-plots-using-ggplot/
qq_plot <- function(plotdata, ldsc=NULL){
    qqplot <- ggplot(plotdata, aes(x = expected, y = observed)) +
        geom_ribbon(aes(ymax = cupper, ymin = clower), fill = "red", alpha = 0.5) +
        geom_point() +
        labs(x = expression(paste("Expected -log"[10],"(", plain(P),")")),
            y = expression(paste("Observed -log"[10],"(", plain(P),")"))) +
        theme_minimal() +
        scale_x_continuous(labels = scales::comma, breaks = scales::pretty_breaks(n = 3)) +
        scale_y_continuous(labels = scales::comma, breaks = scales::pretty_breaks(n = 2)) +
        facet_wrap(~outcome, scales="free_y") +
        theme(
            strip.text.x = element_text(size = 40),
            strip.text = element_text(size = 40),
            axis.text = element_text(size = 40),
            axis.text.x = element_text(size = 40),
            axis.text.y = element_text(size = 40),
            axis.title = element_text(size = 40)
        )
    if(!is.null(ldsc)){
        qqplot <- qqplot +
          geom_text(data = ldsc, aes(label = lab, x = -Inf, y = Inf), hjust = 0, vjust = 1, size=10, parse=T)
    }
    return(qqplot)
}

# load LDSC results
ldsc <- fread("data/ldsc.txt")
names(ldsc) <- c("trait", "intercept", "se")
ldsc$lci <- ldsc$intercept - (1.96 * ldsc$se)
ldsc$uci <- ldsc$intercept + (1.96 * ldsc$se)
ldsc$lab <- paste0("lambda == ", sprintf('%.2f',ldsc$intercept), "  (", sprintf('%.2f',ldsc$lci), "-", sprintf('%.2f',ldsc$uci), ")")
ldsc <- ldsc[ldsc$trait %in% biomarkers]
ldsc$outcome <- sapply(ldsc$trait, function(x) biomarkers_abr[x == biomarkers])

qq_var <- data.frame()
qq_mu <- data.frame()
man_var <- data.frame()
man_mu <- data.frame()
for (i in 1:length(biomarkers)){
    trait_name <- get_trait_name(biomarkers[i])

    # load vGWAS
    data <- get_variants(biomarkers[i])

    # calculate params for QQ plot
    qq_var <- rbind(qq_var, qq_plot_dat(data, "phi_p", biomarkers_abr[i]))
    qq_mu <- rbind(qq_mu, qq_plot_dat(data, "p", biomarkers_abr[i]))
    man_var <- rbind(man_var, data %>% select(chr, pos, phi_p) %>% rename(p=phi_p) %>% mutate(trait = biomarkers_abr[i]))
    man_mu <- rbind(man_var, data %>% select(chr, pos, p) %>% mutate(trait = biomarkers_abr[i]))
}

# qqplot
png("data/gwas_qq_var.log.png", width = 480 * 6, height = 480 * 5)
qq_plot(qq_var)
dev.off()

png("data/gwas_qq_mu.log.png", width = 480 * 6, height = 480 * 5)
qq_plot(qq_mu, ldsc)
dev.off()

png("data/gwas_man_var.log.png", width = 480 * 6, height = 480 * 5)
man_plot(man_var)
dev.off()

png("data/gwas_man_mu.log.png", width = 480 * 6, height = 480 * 5)
man_plot(man_mu)
dev.off()