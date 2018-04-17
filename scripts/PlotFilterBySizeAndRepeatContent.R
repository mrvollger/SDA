#!/usr/bin/env Rscript
library(ggplot2)
suppressPackageStartupMessages(library("argparse"))
# create parser object
parser <- ArgumentParser()
#parser$add_argument("-b", "--bed", default="~/Desktop/data/genomeWide/Mitchell_CHM13_V2/unfiltered.collapses.bed",help="Input bed file")
parser$add_argument("-b", "--bed", default="~/Desktop/data/genomeWide/Mitchell_CHM1/unfiltered.collapses.bed",help="Input bed file")
parser$add_argument("-o", "--png", default="SizeRepeatFilter.png", help="output png file")
parser$add_argument("-s", "--size", type="double", default=9000,help="minimum size of collapse, e.g. 3000")
parser$add_argument("-r", "--repeatContent", type="double", default=75, help="maximum repeat content, e.g. 75")
args <- parser$parse_args()

df = read.table(args$bed, header=F)
contigs = unique(as.character(df[[1]])) 
diffs = c()
for(contig in contigs){
  temp = df[df[1]==contig,]
  len = dim(temp)[1]
  if(len > 1){
    tempdiffs = c(temp[2:len, 3] - temp[1:(len-1), 2],NA)
    diffs = c(diffs, tempdiffs)
  }
  else{
    diffs = c(diffs, NA)
  }
}
#ggplot() + geom_histogram(data=df, aes(V6), bins=50 )

df$diffs = diffs
diffs = diffs[!is.na(diffs)] 

df.melt = df
df.melt$variable = rep("All", length(dim(df.melt)[1]))
df.temp = df.melt[df.melt$V7>=args$size & df.melt$V6<=args$repeatContent,]
df.temp$variable = "Filtered"
df.melt = rbind(df.melt,df.temp)
df.melt$color = rep("#000000", dim(df.melt)[1])
df.melt$color[df.melt$V9 < 50000] = "#990000" 
p=ggplot(data=df.melt, aes(log10(V7),V6)) + 
      geom_density2d(  aes(color=color, alpha=..level..)) + 
      geom_point(alpha = 0.25, aes(color=color)) + 
      scale_y_continuous(breaks = seq(0,100,10), limits=c(0, 100)) +
      facet_wrap( ~ variable, ncol = 2) + xlab("Log10 Length (bp)") + ylab("Percent Repeat Content") +
      theme(legend.position = "none") +
      scale_color_manual(values=c("#990000", "#000000"))
p

ggsave(args$png,plot=p, width = 20, height = 20, units = "cm")

sum(df.melt$variable=="All")
sum(df.melt$variable=="Filtered")
#things close to the end
sum(df.temp$V9 < 50000 )
dim(df.temp)
sum(df$V9 < 50000 & df$V6<=args$repeatContent)
dim(df)


