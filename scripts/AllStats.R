#!/usr/bin/env Rscript
library(ggplot2)
library(scales)
library(RColorBrewer)
library(plyr)
library(gridExtra)
library(grid)
library(gtable)
library(ggrepel)
#library(gtable)
#install.packages("svglite")
library(karyoploteR)
library(GenomicRanges)
library(dplyr)
library(knitr)
library(DT)
library(xtable)
script.dir <- dirname(sys.frame(1)$ofile)
setwd(script.dir)
source("ggplot_theme.R")
bedform = c("chr", "start", "end")


Y= "Yoruban_feb_2018"
C1 = "Mitchell_CHM1_V2"
C13 = "CHM13"
genomes = c(C1, C13, Y)
files = sprintf("~/Desktop/data/genomeWide/%s/LocalAssemblies/localAssemblyStats.tsv", genomes) 
files

df1 =  fread(files[[1]], sep = "\t")
df1$genome = C1
df2 =  fread(files[[2]], sep = "\t")
df2$genome = C13
df3 =  fread(files[[3]], sep = "\t")
df3$genome = Y
df = rbind(df1, df2)
df = rbind(df, df3)

asm = df[df$Status=="Resolved", c("bestChr","bestStart", "bestEnd", "genome")]
names(asm) = c(bedform, "genome")
div = df[df$Status=="Diverged", c("bestChr","bestStart", "bestEnd", "genome")]
names(div) = c(bedform, "genome")
duped = duplicated(df, by=c("collapse") )
dfu = df[ !duped, ]
dfs = rbind(asm, div)
dim(df); dim(asm); dim(div); dim(dfu)


stats <- function(genomex){
  print(genomex)
  asm = asm[asm$genome == genomex,]
  div = div[div$genome == genomex,]
  df = df[df$genome == genomex,]
  dfu = dfu[dfu$genome == genomex,]
  dfs = dfs[dfs$genome == genomex,]
  
  
  
  print(summary(asm$end-asm$start))
  print(summary(dfs$end-dfs$start))
  print(summary(div$end-div$start))
  
  sumres = sum(df[df$Status=="Resolved"]$Length)/10^6
  sumdiv = sum(df[df$Status=="Diverged"]$Length)/10^6
  print(sumres)
  print(sumdiv)
  print(sumres+sumdiv)
  print(sum(dfu$collapseLen)/10^6)
  
  
  asm = reduce( makeGRangesFromDataFrame(asm, ignore.strand = T, keep.extra.columns = T))
  dfs = reduce( makeGRangesFromDataFrame(dfs, ignore.strand = T, keep.extra.columns = T))
  div = reduce( makeGRangesFromDataFrame(div, ignore.strand = T, keep.extra.columns = T))
  
  print( sum(width(asm))/10^6  ) 
  print( sum(width(dfs))/10^6  ) 

}

stats(C1)

stats(C13)

stats(Y)






