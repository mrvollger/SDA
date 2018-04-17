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
suppressPackageStartupMessages(library("argparse"))


library(dplyr)
library(knitr)
library(DT)
library(xtable)


script.dir <- dirname(sys.frame(1)$ofile)
setwd(script.dir)
source("ggplot_theme.R")


readin <- function(file){
    df = read.table(file, header=F)
    colnames(df) = c("chr", "start", "end", "id", "perID", "gene")
    notLOC = grep("LOC|LINC", df[[6]], invert=T)
    df = df[notLOC,]
    
    df$Length = df$end-df$start

    df$Status = "Diverged"
    df$Status[df$perID>99.8] = "Assembled"
    df = unique(df)
    
    return(df)
}


files = Sys.glob("~/Desktop/data/genomeWide/*/LocalAssemblies/betterBlasrMap/all.genes.bed") 
dfs = list()
i = 1
df = NULL
for( file in files ){
  dfs[[i]] = readin(file)
  print(dim(dfs[[i]]))
  df = rbind(df, dfs[[i]])
  i = i + 1
}
dim(df)

acrossAll = Reduce(intersect, list(dfs[[1]]$gene,dfs[[2]]$gene,dfs[[3]]$gene))









makePlot<- function(dfi, dopdf=T, file="~/Desktop/GenicIdeogram.pdf"){
  
  
  # only keep one copy of each gene being mentioned 
  df = dfi[ !duplicated(dfi$gene), ]
  df$chr = as.character(df$chr)
  df$x = (df$start + df$end)/2
  df$y = 0.25 
  # color by name
  df <- df[order(df$gene),]
  df$color = rainbow(length(df$gene))
  df <- df[order(df$chr, df$start),]
  
  assembled = df[df$Status=="Assembled",]
  diverged = df[df$Status=="Diverged",]
  
  tsize = .3
  dopdf=T
  if(dopdf){
    pdf(file, width = 10, height = 30)
  }
  
  kp <- plotKaryotype(  plot.type = 2, genome="hg38")#, chromosomes = c("chr1", "chr2", "chr9"),)
  kpPlotMarkers(kp, chr=assembled$chr, x = assembled$x, y=assembled$y, 
                labels = assembled$gene, line.color=assembled$color, label.color=assembled$color, cex = tsize)
  kpPlotMarkers(kp,  chr=diverged$chr, x = diverged$x, y=diverged$y, 
                labels = diverged$gene, line.color=diverged$color, label.color=diverged$color, data.panel = 2, cex = tsize)
  
  
  
  
  tsize = 0.4
  tt <- gridExtra::ttheme_default(
    core = list(fg_params=list(cex = tsize, parse=T)),
    colhead = list(fg_params=list(cex = tsize, parse=T)))
  

  dfi <- dfi[order(dfi$gene, dfi$perID),]
  rownames(dfi) <- 1:nrow(dfi)
  x = round( length(dfi$gene)/3 )
  df1 = dfi[1:x, c("gene", "perID", "Status")]
  df2 = dfi[x:(2*x), c("gene", "perID", "Status") ]
  df3 = dfi[(2*x):length(dfi$gene), c("gene", "perID", "Status") ]
  print(df1)
  print(df2)
  print(df3)
  #grid.newpage( )
  grid.arrange(tableGrob(df1, theme=tt, rows=NULL), 
               tableGrob(df2, theme=tt, rows=NULL),
               tableGrob(df3, theme=tt, rows=NULL),
               ncol=3)
  #grid.draw(sidebyside)
  #grid.table(df[,c("gene", "perID", "Status")], ) + theme_minimal()
  if(dopdf){
    dev.off()
  }
  
}

makePlot(dfs[[1]], file="~/Desktop/CHM13GenicIdeogram.pdf")

makePlot(dfs[[2]], file="~/Desktop/CHM1GenicIdeogram.pdf")

makePlot(dfs[[3]], file="~/Desktop/YorubanGenicIdeogram.pdf")


for(x in dfs){
  print(dim(x))
}






























file1 = "~/Desktop/data/genomeWide/Mitchell_CHM1_V2/coverage/all.merged.bed"
file2= "~/Desktop/data/genomeWide/Yoruban_feb_2018/coverage/all.merged.bed"
file3= "~/Desktop/data/genomeWide/CHM13/coverage/all.merged.bed"

bed1 = data.table::fread(file1)
bed2 = data.table::fread(file2)
bed3 = data.table::fread(file3)

bed = bed1
colnames(bed) = c("chr", "start", "end", "cov")
covs = sort(bed$cov)
top = covs[ round( .99 * length(covs) ) ]
bot = covs[ round( .01 * length(covs) ) ]
bot
top
withoutcrazy = bed[bed$cov < top & bed$cov > bot, ]
sd(withoutcrazy$cov)
summary(withoutcrazy)


ggplot(withoutcrazy) + geom_histogram(aes(x=cov, y=..density..), position="identity") + 
  geom_density(aes(x=cov, y=..density..))



contig = bed[bed$chr == "LJII01000107.1", ]
contig = contig[contig$start > 5030000, ]

ggplot(contig) + geom_point( aes(start, cov) )

mcov = withoutcrazy$cov
ks.test(mcov,y=rnorm(length(mcov), sd=sd(mcov), mean=mean(mcov) ))
ks.test(mcov,y=runif(length(mcov)))

ks.test(runif(1000),y='pnorm',alternative='two.sided')



