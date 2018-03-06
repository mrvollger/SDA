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


genome = "Mitchell_CHM1_V2"

refdir = "~/Desktop/work/assemblies/hg38/"
asmdir = sprintf("~/Desktop/data/genomeWide/%s/LocalAssemblies/", genome) 
segdir = "~/Desktop/work/assemblies/CHM1/GCA_001297185.1_PacBioCHM1_r2_GenBank_08312015/Segdups/"
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
    pdf(file, width = 15, height = 10)
  }
  if(all){
    kp <- plotKaryotype(genome="hg38", plot.type = 2)
  }else{
    kp <- plotKaryotype(genome="hg38", plot.type = 2, chr = chrs)
  }
  # these makes these tracks easeier to vizualise
  highid <- reduce(highid)
  asm.short <- reduce(asm.short)
  
  # add labels
  kpText(kp, data.panel=1, chr="chr1", x = rep(-15000000,2), y = c(0.25, 0.75), labels = c("Resolved ABP", "Diverged ABP"))
  kpText(kp, data.panel=2, chr="chr1", x = rep(-15000000,2), y = c(0.25,0.75), labels = c("SegDups", "CHM1 Asm"))
  
  # data panael 1, ABP sequences
  kpPlotRegions(kp, data=highid, col=highid.color, layer.margin = 0.05, r0=0, r1=0.45)
  kpPlotRegions(kp, data=lowid, col=lowid.color, layer.margin = 0.05,  r0=0.5, r1=1)
  #kpPlotRegions(kp, data=alreadyDone, col=alreadyDone.color, layer.margin = 0.05, r0=-0.05, r1=0)
  
  # data panel 2, seg dups, colored by status
  kpPlotRegions(kp, data=res, col=res.color, layer.margin = 0.05,  r0=0, r1=0.20, data.panel = 2)
  kpPlotRegions(kp, data=unres, col=unres.color, layer.margin = 0.05, r0=0.25, r1=.45, data.panel = 2)
  
  # data panel 2, CHM1 assembly plot
  #asm.short <- reduce(resolvedSegdups)
  #asm.long <- reduce(unresolvedSegdups)
  kpPlotRegions(kp, data=asm.long, col=asm.long.color, layer.margin = 0.05, r0=.70, r1=0.5, data.panel = 2)
  kpPlotRegions(kp, data=asm.short, col=asm.short.color, layer.margin = 0.05, r0=.95, r1=0.75, data.panel = 2)
  
  if(dopdf){
    dev.off()
  }
}

#
# plot commands 
#
subset = c("chr1","chr2", "chr9", "chr16", "chr14")
plotKP(chrs=subset)

if(T){
  subset = c("chr1","chr2", "chr9", "chr16", "chr14")
  plotKP(subset, dopdf=T, name="subset")
  plotKP(all=T, dopdf=T, name="ideogram")
  for(chr in levels(seqnames(asm))){
    plotKP(chr, dopdf=T, name=chr)
  }
}























if(F){
#
# not used 
#

#
# intersection analysis
#

# tests
#asm <- IRanges(c(3), c(10))
#segdups <- IRanges(c(4, 1), c(7, 12))
asmNotInSegDup = subsetByOverlaps(asm, segdups, type=c("within"), invert=T)
asmInSegDup = subsetByOverlaps(asm, segdups, type=c("within"), invert=F)

resolvedSegdups = subsetByOverlaps(segdups, asmNotInSegDup, type="within")
unresolvedSegdups = subsetByOverlaps(segdups, asmNotInSegDup, type="within", invert=T)

length(asm); length(asmNotInSegDup); length(asmInSegDup); length(resolvedSegdups); length(unresolvedSegdups)
sum( width(reduce(unresolvedSegdups)) )/10^6
sum( width(reduce(resolvedSegdups)) )/10^6

# require not being part of a larger seg dup
sum(width(reduce(highid)))/10^6
abpNotInSegDup = subsetByOverlaps(highid, unresolvedSegdups, type=c("within"), invert=T)
sum(width(reduce(abpNotInSegDup)))/10^6
segDupsInABP = subsetByOverlaps(unresolvedSegdups, abpNotInSegDup, type=c("within"))
sum(width(reduce(segDupsInABP)))/10^6
# again without the requiremnt 
segDupsInABP = subsetByOverlaps(unresolvedSegdups, (highid), type=c("within"))
sum(width(reduce(segDupsInABP)))/10^6



rescaleDecimal = function(points, botDecimal, topDecimal){
  points[]
  ints = floor(points)
  decimals = points - ints
  scaled = decimals * (topDecimal - botDecimal)/(1 - 0) + botDecimal
  return(ints + scaled)
}

status = c(rep("Resolved", length(resolvedSegdups)), 
           rep("Unresolved", length(unresolvedSegdups)),
           rep("New", length(segDupsInABP)))
x = c(log10(width(resolvedSegdups)), log10(width(unresolvedSegdups)), log10(width(segDupsInABP)))
y = c(resolvedSegdups$perID, unresolvedSegdups$perID, segDupsInABP$perID)*100
tempdf = data.frame(x,y,status )
tempdf$notLog = 10^tempdf$x
tempdf[tempdf$y <= 90.0, "y"] = 90.00001
tempdf[tempdf$y >= 100.00, "y"] = 99.99999
tempdf$bin = cut( tempdf$y, breaks = seq(90,100))

unique(tempdf$bin)

bybin = split(tempdf, tempdf$bin)
for(bin in bybin){
  R = sum( bin[bin$status == "Resolved",]$notLog  )
  U = sum( bin[bin$status == "Unresolved",]$notLog  )
  N = sum( bin[bin$status == "New",]$notLog )
  print(R/(R+U) )
  print( (R+N)/(R+U) )
  print(as.character( bin$bin[1]))
}

tempdf$dec = c(  rescaleDecimal( tempdf$y[tempdf$status=="Resolved"], 0.0, .33), 
                 rescaleDecimal( tempdf$y[tempdf$status=="Unresolved"], 0.33, .66),
                 rescaleDecimal( tempdf$y[tempdf$status=="New"], 0.66, 1))
ggplot(tempdf)+geom_point(aes(x=x, y=dec, color = status), alpha = .5) 

ggplot(tempdf) +
  geom_density_2d(aes(x=x, y=y, color=status, alpha=(..level..) ),
                  size = 2, n = 50, binwidth=.015) +
  theme_classic()


print("density")

}

