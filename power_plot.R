library("ggplot2")
library("dplyr")
library("tidyr")
library("broom")
library("data.table")

#' Function to plot the power analysis
#' @param h1 dataframe from power_funs.R:calc_power
#' @param title Title of plot
#' @param x_title X axis title
#' @param x_name Name of the X variable
#' @param y_name Name of the Y variable
#' @param ymin_name Name of the Y variable lower 95% CI
#' @param ymax_name Name of the Y variable upper 95% CI
#' @param xh Rotate the X axis text 90o
#' @param group Name of the series
#' @param color Name of the variable to color the series
plot_power <- function(h1, title, x_title, x_name, y_name, ymin_name, ymax_name, xh=FALSE, group=1, color=NULL, alpha=0.05){
    p <- ggplot(data=h1, aes_string(x=x_name, y=y_name, ymin=ymin_name, ymax=ymax_name, group=group, color=color)) +
        geom_line() + 
        geom_point() + 
        geom_errorbar(width=.1, position=position_dodge(0.2)) +
        theme_classic() + 
        ggtitle(title) +
        xlab(x_title) + 
        ylab(paste0("Power (alpha=", alpha, ")")) +
        scale_y_continuous(limits = c(0, 1), breaks = scales::pretty_breaks(n = 10))
    
    if (xh){
        p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    }

    return(p)
}

#' Function to estimate the power of an MC experiment
#' @param results Dataframe containing: P value and analysis group(s)
#' @param field Name of P value to summarise
#' @param n_sim Number of simulations performed
#' @param grp_name A vector of fields for use in grouping the analysis
#' @param alpha P threshold
calc_power <- function(results, field, n_sim, grp_name, alpha=0.05){
    # threshold P value
    results <- results %>%
        mutate(pos = as.numeric(get(field) < alpha))
    
    # count positives and estimate power with 95% CI
    h1 <- results %>%
        select(pos, all_of(grp_name)) %>%
        drop_na() %>%
        group_by_at(vars(all_of(grp_name))) %>%
        summarise(h1 = sum(pos)) %>%
        rowwise() %>%
        mutate(est_power = tidy(binom.test(h1, n_sim))$estimate) %>% 
        mutate(est_power_low = tidy(binom.test(h1, n_sim))$conf.low) %>%
        mutate(est_power_high = tidy(binom.test(h1, n_sim))$conf.high)

    # factorise grouping variables
    h1[grp_name] <- lapply(h1[grp_name], factor)

    return(h1)
}


# load data
d <- fread("results.csv")

# process data
p <- calc_power(d, "P.cpp", 200, c("phi", "lambda"))

# plot
pdf("plot.pdf")
plot_power(p, "vGWAS power to detect GxE effect", "Sample size inflation factor", "lambda", "est_power", "est_power_low", "est_power_high", group="phi", color="phi")
dev.off()