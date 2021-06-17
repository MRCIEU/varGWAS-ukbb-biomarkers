library('data.table')
library('dplyr')
library('stringr')
library('ggplot2')
source("funs.R")
set.seed(1234)

data <- data.frame()
for (i in 1:length(biomarkers)){
  # extract trait name
  trait_name <- get_trait_name(biomarkers[i])

  # load phenotype
  d <- fread(paste0("data/", biomarkers[i], ".txt"), select=biomarkers[i], col.names=c("y"))

  # add trait name & append
  d$trait <- biomarkers_abr[i]
  data <- rbind(data, d)
}

# histogram
png("data/hist.png", width = 480 * 1.5, height = 480 * 1.25)
ggplot(data, aes(x=y)) +
  geom_histogram() +
  theme_minimal() +
  scale_x_continuous(labels = scales::comma, breaks = scales::pretty_breaks(n = 3)) +
  scale_y_continuous(labels = scales::comma, breaks = scales::pretty_breaks(n = 3)) +
  facet_wrap(~trait, scales="free") +
  theme(
    axis.title=element_blank()
)
dev.off()

# qqplot
png("data/qq.png", width = 480 * 1.5, height = 480 * 1.25)
ggplot(data, aes(sample=y)) +
  stat_qq() +
  stat_qq_line() +
  theme_minimal() +
  scale_x_continuous(labels = scales::comma, breaks = scales::pretty_breaks(n = 3)) +
  scale_y_continuous(labels = scales::comma, breaks = scales::pretty_breaks(n = 3)) +
  facet_wrap(~trait, scales="free") +
  theme(
    axis.title=element_blank()
)
dev.off()