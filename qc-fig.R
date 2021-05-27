library('data.table')
library('dplyr')
library('stringr')
library('ggpubr')
source("funs.R")
set.seed(1234)

mean_man <- list()
mean_qq <- list()
var_man <- list()
var_qq <- list()
labs <- rep(NA, length(biomarkers))

# Taken from https://danielroelfs.com/blog/how-i-create-manhattan-plots-using-ggplot/
man_plot <- function(gwas_data, pval_col, sig=5e-8/30, ylim=30){
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

    gwas_data$lp <- -log10(gwas_data[[pval_col]])

    manhplot <- ggplot(gwas_data, aes(x = bp_cum, y = lp, color = as_factor(chr), size = lp)) +
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
            panel.border = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()
        )

    return(manhplot)
}

# Taken from https://danielroelfs.com/blog/how-i-make-qq-plots-using-ggplot/
qq_plot <- function(sumstats.data, pval_col, ci=0.95){
    nSNPs <- nrow(sumstats.data)
    plotdata <- data.frame(
        observed = -log10(sort(sumstats.data[[pval_col]])),
        expected = -log10(ppoints(nSNPs)),
        clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = seq(nSNPs), shape2 = rev(seq(nSNPs)))),
        cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = seq(nSNPs), shape2 = rev(seq(nSNPs))))
    )
    qqplot <- ggplot(plotdata, aes(x = expected, y = observed)) +
        geom_ribbon(aes(ymax = cupper, ymin = clower), fill = "grey30", alpha = 0.5) +
        geom_step(color = "#183059", size = 1.1, direction = "vh") +
        geom_segment(data = . %>% filter(expected == max(expected)), 
                    aes(x = 0, xend = expected, y = 0, yend = expected),
                    size = 1.25, alpha = 0.5, color = "grey30", lineend = "round") +
        labs(x = expression(paste("Expected -log"[10],"(", plain(P),")")),
            y = expression(paste("Observed -log"[10],"(", plain(P),")"))) +
        theme_minimal() +
        theme(
            axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()
        )
    return(qqplot)
}

for (i in 1:length(biomarkers)){
    trait_name <- biomarkers[i]
    trait_name <- str_split(trait_name, "\\.", simplify = TRUE)[,1]
    trait_name <- gsub("_", " ", trait_name)
    trait_name <- str_to_title(trait_name)
    labs[i] <- trait_name
    message(paste0("trait ", biomarkers[i]))
    message(paste0("trait name ", trait_name))

    # load vGWAS and QC
    data <- get_variants(biomarkers[i])

    # manhattan
    mean_man[[i]] <- man_plot(data, "p")
    var_man[[i]] <- man_plot(data, "phi_p")

    # qq plot
    mean_qq[[i]] <- qq_plot(data, "p")
    var_qq[[i]] <- qq_plot(data, "phi_p")
}

png("var_qq.png", width = 480 * 2.5, height = 480 * 3)
p <- ggarrange(plotlist=var_qq, labels = biomarkers_abr, ncol = 5, nrow = 6, align = "hv", font.label=list(size = 18))
annotate_figure(p, left = grid::textGrob("Observed -log10(P)", rot = 90, vjust = 1, gp = grid::gpar(cex = 1.5)),
                    bottom = grid::textGrob("Expected -log10(P)", gp = grid::gpar(cex = 1.5)))
dev.off()

png("mean_qq.png", width = 480 * 2.5, height = 480 * 3)
p <- ggarrange(plotlist=mean_qq, labels = biomarkers_abr, ncol = 5, nrow = 6, align = "hv", font.label=list(size = 18))
annotate_figure(p, left = grid::textGrob("Observed -log10(P)", rot = 90, vjust = 1, gp = grid::gpar(cex = 1.5)),
                    bottom = grid::textGrob("Expected -log10(P)", gp = grid::gpar(cex = 1.5)))
dev.off()

png("mean_man.png", width = 480 * 2.5, height = 480 * 3)
p <- ggarrange(plotlist=mean_man, labels = biomarkers_abr, ncol = 5, nrow = 6, align = "hv", font.label=list(size = 18))
annotate_figure(p, left = grid::textGrob("-log10(P)", rot = 90, vjust = 1, gp = grid::gpar(cex = 1.5)),
                    bottom = grid::textGrob("Chromosome", gp = grid::gpar(cex = 1.5)))
dev.off()

png("var_man.png", width = 480 * 2.5, height = 480 * 3)
p <- ggarrange(plotlist=var_man, labels = biomarkers_abr, ncol = 5, nrow = 6, align = "hv", font.label=list(size = 18))
annotate_figure(p, left = grid::textGrob("-log10(P)", rot = 90, vjust = 1, gp = grid::gpar(cex = 1.5)),
                    bottom = grid::textGrob("Chromosome", gp = grid::gpar(cex = 1.5)))
dev.off()