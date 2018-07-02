#!/usr/bin/env Rscript

library(ggplot2)
library(scales)
library(RColorBrewer)
library(plyr)
library(gridExtra)
library(ggrepel)
library(grid)
library(gtable)
#source("http://bioconductor.org/biocLite.R")
#biocLite("karyoploteR")
library(karyoploteR)
library(GenomicRanges)
suppressPackageStartupMessages(library("argparse"))
bedform = c("chr", "start", "end")

genome = "CHM13"
genome = "Mitchell_CHM1_V2"
genome = "Mitchell_CHM1"
genome = "Yoruban_feb_2018"


refdir = "~/Desktop/work/assemblies/hg38/"
asmdir = sprintf("~/Desktop/data/genomeWide/%s/LocalAssemblies/", genome) 
segdir <- Sys.glob(sprintf("~/Desktop/work/assemblies/%s/*/Segdups/", genome) )[1]
plotsdir = sprintf("~/Desktop/data/genomeWide/%s/plots/", genome) 

# create parser object
parser <- ArgumentParser()
parser$add_argument("-r", "--refdir", default=refdir, help="Input tsv file")
parser$add_argument("-a", "--asmdir", default=asmdir, help="destination folder")
parser$add_argument("-s", "--segdir", default=segdir, help="list of resolved seg dups")
parser$add_argument("-p", "--plotsdir", default=plotsdir, help=" list of unresolved seg dups")
args <- parser$parse_args()

refdir = args$refdir
asmdir = args$asmdir 
segdir = args$segdir
plotsdir = args$plotsdir


if( ! dir.exists(plotsdir)){
  dir.create(plotsdir)
}
if(! dir.exists(paste0(plotsdir, "ideogram")) ){
  dir.create(paste0(plotsdir, "ideogram"))
}


# read in all segdups, unmerged
file = paste0(refdir, "ucsc.unmerged.segdups.bed")
segdups = read.table(file, header=F)
names(segdups) = c(bedform, "perID")
segdups = makeGRangesFromDataFrame(segdups, ignore.strand = T, keep.extra.columns = T)
segdups.color = "#000000"

# read in high ID sequences 
file = paste0(asmdir, "betterBlasrMap/highID.bed" )
highid = read.table(file, header=F)
names(highid) = c(bedform, "perID")
highid = makeGRangesFromDataFrame(highid, ignore.strand = T, keep.extra.columns = T)
highid.color = "#000000"

# read in low ID sequences 
file = paste0(asmdir, "betterBlasrMap/lowID.bed" )
lowid = read.table(file, header=F)
names(lowid) = c(bedform, "perID")
lowid = makeGRangesFromDataFrame(lowid, ignore.strand = T, keep.extra.columns = T)
lowid.color = "#A9A9A9"

# read in resolved seg dups
file = paste0(segdir, "asm.mean.resolved")
res = read.table(file, header=F)
names(res) = c(bedform, "perID")
res = makeGRangesFromDataFrame(res, ignore.strand = T, keep.extra.columns = T)
res.color = "#000000"

# read in unresolved seg dups
file = paste0(segdir, "asm.mean.unresolved")
unres = read.table(file, header=F)
names(unres) = c(bedform, "perID")
unres = makeGRangesFromDataFrame(unres, ignore.strand = T, keep.extra.columns = T)
unres.color = "#b20000"


# read in unresolved seg dups
file = paste0(asmdir, "all.ref.fasta.bed")
collapsed = read.table(file, header=F)
names(collapsed) = c(bedform)
collapsed = makeGRangesFromDataFrame(collapsed, ignore.strand = T, keep.extra.columns = T)
collapsed <- reduce(collapsed)
collapsed.color = "#b20000"


# read in assembly 
file = paste0(segdir, "asm.bed5")
asm = read.table(file, header=F)
names(asm) = c(bedform, "conitg", "mapq")
cond = (asm$end - asm$start) > 1000000
asm.long = asm[cond, ]
asm.long = makeGRangesFromDataFrame(asm.long, ignore.strand = T, keep.extra.columns = T)
asm.long.color = "#0000A0"
asm.short = asm[!cond,]
asm.short = makeGRangesFromDataFrame(asm.short, ignore.strand = T, keep.extra.columns = T)
asm.short.color = "#ADD8E6"
asm = makeGRangesFromDataFrame(asm, ignore.strand = T, keep.extra.columns = T)
asm.color = "#000000"

# intersect of short asm and unresolved
alreadyDone = GenomicRanges::intersect(unres, asm.short, ignore.strand=T)
alreadyDone.color = "#ADD8E6"




#
#  the plot function
#
plotKP <- function(chrs, dopdf=FALSE, name="ideogram", all=FALSE){
  if(dopdf){
    file = paste0(paste0(plotsdir,"ideogram/"), name)
    file = paste0(file, ".pdf")
    pdf(file, width = 10, height = 10 )
  }
  if(all){
    kp <- plotKaryotype(genome="hg38", plot.type = 2)
  }else{
    kp <- plotKaryotype(genome="hg38", plot.type = 2, chr = chrs)
  }
  # these makes these tracks easeier to vizualise
  lowid <- reduce(lowid)
  highid <- reduce(highid)
  asm.short <- reduce(asm.short)
  
  # add labels
  #kpText(kp, data.panel=1, chr="chr1", x = rep(-15000000,2), y = c(0.25, 0.75), labels = c("Assembled ABP", "Diverged ABP"))
  #kpText(kp, data.panel=2, chr="chr1", x = rep(-15000000,2), y = c(0.10, 0.35, 0.75), labels = c("Collapsed","SegDups", "DeNovo Asm"))
  
  # data panael 1, ABP sequences
  #kpPlotRegions(kp, data=highid, col=highid.color, layer.margin = 0.05, r0=0, r1=0.45)
  #kpPlotRegions(kp, data=lowid, col=lowid.color, layer.margin = 0.05,  r0=0.5, r1=1)
  kpPlotRegions(kp, data=lowid, col=asm.short.color, layer.margin = 0.05,  r0=0.5-.5, r1=1)
  kpPlotRegions(kp, data=highid, col=asm.long.color, layer.margin = 0.05, r0=0, r1=0.45+.55)
  # data panel 2, seg dups, colored by status
  #kpPlotRegions(kp, data=res, col=res.color, layer.margin = 0.05, r0=0.25, r1=.45, data.panel = 2)
  #kpPlotRegions(kp, data=unres, col=unres.color, layer.margin = 0.05, r0=0.25, r1=.45, data.panel = 2)
  #kpPlotRegions(kp, data=collapsed, col=collapsed.color, layer.margin = 0.05,  r0=0, r1=0.20, data.panel = 2)
  
  kpPlotRegions(kp, data=res, col=res.color, layer.margin = 0.05, r0=0.5-.5, r1=1, data.panel = 2)
  kpPlotRegions(kp, data=unres, col=unres.color, layer.margin = 0.05, r0=0.5-.5, r1=1, data.panel = 2)
  #kpPlotRegions(kp, data=collapsed, col=lowid.color, layer.margin = 0.05,  r0=0, r1=0.45, data.panel = 2)
  
  
  # data panel 2, CHM1 assembly plot
  #kpPlotRegions(kp, data=asm.long, col=asm.long.color, layer.margin = 0.05, r0=.70, r1=0.5, data.panel = 2)
  #kpPlotRegions(kp, data=asm.short, col=asm.short.color, layer.margin = 0.05, r0=.95, r1=0.75, data.panel = 2)
  
  if(dopdf){
    dev.off()
  }
}

#
# plot commands 
#
subset = c("chr1", "chr16")
plotKP(chrs=subset)
if(T){
  subset = c("chr1", "chr16")
  plotKP(subset, dopdf=T, name="subset")
  plotKP(all=T, dopdf=T, name="ideogram")
  for(chr in levels(seqnames(asm))){
    #plotKP(chr, dopdf=T, name=chr)
    print("x")
  }
}

if(T){
sds = c(res, unres); sds = sds[ width(sds) > 50000,]
kp <- plotKaryotype(genome="hg38", plot.type = 2, chr = chrs)
kpPlotRegions(kp, data=asm.long, col=asm.long.color, layer.margin = 0.05, r0=0, r1=.33, data.panel = 2)
kpPlotRegions(kp, data=reduce(asm.short), col=asm.short.color, layer.margin = 0.05, r0=.33, r1=0.66, data.panel = 2)

kpPlotRegions(kp, data=sds, col=unres.color, layer.margin = 0.05, r0=.66, r1=1, data.panel = 2)


}

sum(width(reduce(highid)))/10^6
sum(width(reduce(lowid)))/10^6
sum(width(reduce(c(lowid, highid))))/10^6
summary(width(highid))

summary(width(lowid))



#
# combined plot
#
if(F){

x = reduce( c( highid, lowid  ) )
sum(width(x))/10^6

genome = "Yoruban_feb_2018"
asmdir = sprintf("~/Desktop/data/genomeWide/%s/LocalAssemblies/", genome) 
# read in high ID sequences 
file = paste0(asmdir, "betterBlasrMap/highID.bed" )
highid = read.table(file, header=F)
names(highid) = c(bedform, "perID")
yhi =  reduce(makeGRangesFromDataFrame(highid, ignore.strand = T, keep.extra.columns = T))
# read in low ID sequences 
file = paste0(asmdir, "betterBlasrMap/lowID.bed" )
lowid = read.table(file, header=F)
names(lowid) = c(bedform, "perID")
ylow = reduce( makeGRangesFromDataFrame(lowid, ignore.strand = T, keep.extra.columns = T))


genome = "CHM13"
asmdir = sprintf("~/Desktop/data/genomeWide/%s/LocalAssemblies/", genome) 
# read in high ID sequences 
file = paste0(asmdir, "betterBlasrMap/highID.bed" )
highid = read.table(file, header=F)
names(highid) = c(bedform, "perID")
chm13hi =  reduce(makeGRangesFromDataFrame(highid, ignore.strand = T, keep.extra.columns = T))
# read in low ID sequences 
file = paste0(asmdir, "betterBlasrMap/lowID.bed" )
lowid = read.table(file, header=F)
names(lowid) = c(bedform, "perID")
chm13low = reduce( makeGRangesFromDataFrame(lowid, ignore.strand = T, keep.extra.columns = T))

genome = "Mitchell_CHM1_V2"
asmdir = sprintf("~/Desktop/data/genomeWide/%s/LocalAssemblies/", genome) 
# read in high ID sequences 
file = paste0(asmdir, "betterBlasrMap/highID.bed" )
highid = read.table(file, header=F)
names(highid) = c(bedform, "perID")
chm1hi =  reduce(makeGRangesFromDataFrame(highid, ignore.strand = T, keep.extra.columns = T))
# read in low ID sequences 
file = paste0(asmdir, "betterBlasrMap/lowID.bed" )
lowid = read.table(file, header=F)
names(lowid) = c(bedform, "perID")
chm1low = reduce( makeGRangesFromDataFrame(lowid, ignore.strand = T, keep.extra.columns = T))

  
  
  
    
lowid.color = "#A9A9A9"
highid.color = "#000000"



offset = -0.08
kp <- plotKaryotype(genome="hg38", plot.type = 1, chr = subset)
# data panael 1, ABP sequences
kpPlotRegions(kp, data=chm1hi, col=highid.color, layer.margin = 0.05, r0=0.00 +offset, r1=0.15 +offset)
kpPlotRegions(kp, data=chm13hi, col=highid.color, layer.margin = 0.05,  r0=0.18 +offset, r1=0.33 +offset)
kpPlotRegions(kp, data=yhi, col=highid.color, layer.margin = 0.05, r0=.36 +offset, r1=0.51 +offset)

kpPlotRegions(kp, data=chm1low, col=lowid.color, layer.margin = 0.05,  r0=.54 +offset, r1=0.69 +offset)
kpPlotRegions(kp, data=chm13low, col=lowid.color, layer.margin = 0.05, r0=.72 +offset, r1=0.87 +offset)
kpPlotRegions(kp, data=ylow, col=lowid.color, layer.margin = 0.05,  r0=.90 +offset, r1=1.05 +offset)


file = paste0("~/Desktop/", "compare.pdf")
pdf(file, width = 15, height = 30)
offset = -0.08
kp <- plotKaryotype(genome="hg38", plot.type = 1)
# data panael 1, ABP sequences
kpPlotRegions(kp, data=chm1hi, col=highid.color, layer.margin = 0.05, r0=0.00 +offset, r1=0.15 +offset)
kpPlotRegions(kp, data=chm13hi, col=highid.color, layer.margin = 0.05,  r0=0.18 +offset, r1=0.33 +offset)
kpPlotRegions(kp, data=yhi, col=highid.color, layer.margin = 0.05, r0=.36 +offset, r1=0.51 +offset)

kpPlotRegions(kp, data=chm1low, col=lowid.color, layer.margin = 0.05,  r0=.54 +offset, r1=0.69 +offset)
kpPlotRegions(kp, data=chm13low, col=lowid.color, layer.margin = 0.05, r0=.72 +offset, r1=0.87 +offset)
kpPlotRegions(kp, data=ylow, col=lowid.color, layer.margin = 0.05,  r0=.90 +offset, r1=1.05 +offset)
dev.off()




offset = -0.08
kp <- plotKaryotype(genome="hg38", plot.type = 1, chr = subset)
# data panael 1, ABP sequences
kpPlotRegions(kp, data=reduce(c(chm1hi,chm1low)), col=highid.color, layer.margin = 0.05, r0=0.00 +offset, r1=0.15 +offset)
kpPlotRegions(kp, data=reduce(c(chm13hi,chm13low)), col=highid.color, layer.margin = 0.05,  r0=0.18 +offset, r1=0.33 +offset)
kpPlotRegions(kp, data=reduce(c(yhi,ylow)), col=highid.color, layer.margin = 0.05, r0=.36 +offset, r1=0.51 +offset)

kpPlotRegions(kp, data=, col=lowid.color, layer.margin = 0.05,  r0=.54 +offset, r1=0.69 +offset)
kpPlotRegions(kp, data=, col=lowid.color, layer.margin = 0.05, r0=.72 +offset, r1=0.87 +offset)
kpPlotRegions(kp, data=, col=lowid.color, layer.margin = 0.05,  r0=.90 +offset, r1=1.05 +offset)


}

