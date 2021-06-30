library('ggplot2')
library("ggExtra")
set.seed(123)
setEPS()

n_obs <- 1000
x <- rbinom(n_obs, 2, .5)
u <- runif(n_obs)
y0 <- 20*x + 0*x*u + rnorm(n_obs)
y1 <- 20*x + 20*x*u + rnorm(n_obs)
df <- data.frame(rbind(
    data.frame(
        x=x,
        y=scale(y0),
        d="A"
    ),
    data.frame(
        x=x,
        y=scale(y1),
        d="B"
    )
))
df$x <- as.factor(df$x)

p1 <- ggplot(df, aes(x=x, y=y)) +
    geom_boxplot() +
    xlab("SNP") +
    ylab("Outcome") +
    theme_classic() +
    theme(
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.background=element_blank(),
        strip.text.x = element_text(size = 20),
        text = element_text(size=20)
    )
p1 <- p1 + facet_grid(.~df$d)

ggsave("infograpic.eps", p1)