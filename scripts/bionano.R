#!/usr/bin/env Rscript
library(ggplot2)
library(scales)
library(RColorBrewer)
library(plyr)
library(gridExtra)
library(ggrepel)
library(grid)
library(gtable)

file = "~/Desktop/lowid.txt"
file2 = "~/Desktop/highid.txt"
lowid = read.table(file, header = T)
lowid$id="low"
highid = read.table(file2, header = T)
highid$id="high"

df = rbind(lowid, highid)
df = df[ df$altWeightedAvgCov > 0.0 & df$hg38FragmentWeightedAvgCov > 0.0 , ]
df$ratio = (  df$altWeightedAvgCov/df$hg38FragmentWeightedAvgCov  )
df$logratio = log10(  df$ratio )

ggplot(df) + geom_boxplot(aes(x=id, y=logratio, color=id)) 

sum(lowid$hg38FragmentWeightedAvgCov < lowid$altWeightedAvgCov)
sum(lowid$hg38FragmentWeightedAvgCov > lowid$altWeightedAvgCov)
sum(highid$hg38FragmentWeightedAvgCov < highid$altWeightedAvgCov)
sum(highid$hg38FragmentWeightedAvgCov > highid$altWeightedAvgCov)


x = df[df$id == "low", ]$ratio
y = df[df$id == "high", ]$ratio

wilcox.test(highid$altWeightedAvgCov, highid$hg38FragmentWeightedAvgCov, paired=TRUE)
wilcox.test(lowid$altWeightedAvgCov, lowid$hg38FragmentWeightedAvgCov, paired=TRUE)


ggplot(df, aes(hg38FragmentWeightedAvgCov, altWeightedAvgCov)) + geom_point(aes(color=id)) + # stat_summary(fun.data=mean_cl_normal) + 
  geom_smooth(method='lm',formula=y~x) +
  coord_fixed() + geom_abline(slope=1, intercept=0, color = "red") + facet_grid(id ~ .)




ggplot(highid, aes(hg38FragmentWeightedAvgCov, altWeightedAvgCov)) + geom_point() + # stat_summary(fun.data=mean_cl_normal) + 
  #geom_smooth(method='lm',formula=y~x) +
  coord_fixed() + geom_abline(slope=1, intercept=0, color = "red") 

ggplot(lowid, aes(hg38FragmentWeightedAvgCov, altWeightedAvgCov)) + geom_point() + # stat_summary(fun.data=mean_cl_normal) + 
  #geom_smooth(method='lm',formula=y~x) +
  coord_fixed()+ geom_abline(slope=1, intercept=0, color = "red")



sum(df$altWeightedAvgCov)
sum(df$hg38FragmentWeightedAvgCov)

sum(lowid$altWeightedAvgCov)
sum(lowid$hg38FragmentWeightedAvgCov)

