setwd("~/projects/OptiMimi-j5/test")

library(reshape2)
library(plyr)

attempts <- read.csv("attempts-compare.csv")
attempts$samplemult <- c("1x", rep(c("1x", "2x", "4x"), each=9))

df <- melt(attempts[,c(1:12, 23:24)], id.vars=c('mcperlife', 'samplemult', 'algorithm', 'time'))

levels(df$variable) <- as.character(1:10)
df$algorithm <- revalue(df$algorithm, c("social"="Homegrown", "biological"="BBO SNES", "sampled"="Sampled"))
df$algorithm <- factor(df$algorithm, c("Homegrown", "BBO SNES", "Sampled"))

df.exact <- rbind(subset(df, is.na(algorithm)), subset(df, is.na(algorithm)), subset(df, is.na(algorithm)))
df.exact$algorithm <- rep(c("Homegrown", "BBO SNES", "Sampled"), each=10)

library(ggplot2)

ggplot(subset(df, algorithm != 'exact'), aes(variable, value, colour=factor(mcperlife), linetype=factor(samplemult), group=paste0(mcperlife, samplemult, algorithm))) +
    facet_grid(algorithm ~ .) +
    #geom_hline(yintercept=1/3) +
    geom_line() + scale_colour_discrete(name="Monte Carlos / Objective: ") +
    geom_line(data=df.exact, colour="black", linetype="solid") +
    scale_linetype_manual(name="Evaluations / Monte Carlo: ", breaks=c('1x', '2x', '4x'), values=c('dotted', 'dashed', 'solid')) +
    scale_y_continuous(limits=c(0, 1), expand=c(0, 0)) +
    theme_bw() + theme(legend.position="bottom") + xlab("Period") +
    ylab("Consumption level") + theme(panel.spacing = unit(1, "lines"))
#ggsave(paste0("uncertainty", suffix, ".pdf"), width=5, height=3.5)

attempts$rmse <- sqrt(((attempts$cons1 - attempts$cons1[1])^2 + (attempts$cons2 - attempts$cons2[1])^2 + (attempts$cons3 - attempts$cons3[1])^2 + (attempts$cons4 - attempts$cons4[1])^2 + (attempts$cons5 - attempts$cons5[1])^2 + (attempts$cons6 - attempts$cons6[1])^2 + (attempts$cons7 - attempts$cons7[1])^2 + (attempts$cons8 - attempts$cons8[1])^2 + (attempts$cons9 - attempts$cons9[1])^2 + (attempts$cons10 - attempts$cons10[1])^2) / 10)

attempts$algorithm <- revalue(attempts$algorithm, c("social"="Homegrown", "biological"="BBO SNES", "sampled"="Sampled"))
attempts$algorithm <- factor(attempts$algorithm, c("Homegrown", "BBO SNES", "Sampled"))

ggplot(subset(attempts, algorithm != 'exact'), aes(time, rmse)) +
    facet_grid(. ~ algorithm) +
    geom_point() + geom_smooth(method='lm') +
    scale_y_log10() + theme_bw() +
    xlab("Computation time (s)") + ylab("RMSE from true optimum")
