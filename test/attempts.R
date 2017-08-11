setwd("~/projects/OptiMimi-j5/test")

library(reshape2)

attempts <- read.csv("attempts-full.csv")

grouped <- attempts[1, ]
for (minmc in seq(1, 40, by=15))
    grouped <- rbind(grouped, colMeans(attempts[attempts$mcperlife %in% (0:14 + minmc),]))
for (minmc in seq(50, 100, by=30))
    grouped <- rbind(grouped, colMeans(attempts[attempts$mcperlife %in% (0:29 + minmc),]))
for (minmc in seq(100, 200, by=40))
    grouped <- rbind(grouped, colMeans(attempts[attempts$mcperlife %in% (0:39 + minmc),]))

df <- melt(grouped, id.vars='mcperlife')

levels(df$variable) <- as.character(1:10)

library(ggplot2)

ggplot(df, aes(variable, value, colour=mcperlife, group=mcperlife)) +
    geom_line() + scale_colour_gradientn(name="Monte Carlos per evaluation: ", colours=rainbow(4)) +
    theme_bw() + theme(legend.position="bottom") + xlab("Period") +
    ylab("Consumption level")
ggsave("uncertainty-full.pdf", width=5, height=3.5)
