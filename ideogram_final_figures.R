library(ggplot2)
library(scales)
library(RColorBrewer)
library(plyr)
library(gridExtra)
library(ggrepel)
library(grid)
library(gtable)
library(karyoploteR)
library(GenomicRanges)
bedform = c("chr", "start", "end")


# read in high ID sequences 
file = "~/Desktop/data/genomeWide/Mitchell_CHM1/LocalAssemblies/betterBlasrMap/highID.bed"
highid = read.table(file, header=F)
names(highid) = c(bedform, "perID")
highid = makeGRangesFromDataFrame(highid, ignore.strand = T, keep.extra.columns = T)
highid.color = "#000000"
highid = reduce(highid)

# read in low ID sequences 
file = "~/Desktop/data/genomeWide/Mitchell_CHM1/LocalAssemblies/betterBlasrMap/lowID.bed"
lowid = read.table(file, header=F)
names(lowid) = c(bedform, "perID")
lowid = makeGRangesFromDataFrame(lowid, ignore.strand = T, keep.extra.columns = T)
lowid.color = "#A9A9A9"
#lowid = reduce(lowid)

# read in resolved seg dups
file = "~/Desktop/data/genomeWide/Mitchell_CHM1/segdups/CHM1_Falcon_P6_27.fasta.resolved"
res = read.table(file, header=F)
names(res) = c(bedform, "perID")
res = makeGRangesFromDataFrame(res, ignore.strand = T, keep.extra.columns = T)
res.color = "#000000"

# read in unresolved seg dups
file = "~/Desktop/data/genomeWide/Mitchell_CHM1/segdups/CHM1_Falcon_P6_27.fasta.unresolved"
unres = read.table(file, header=F)
names(unres) = c(bedform, "perID")
unres = makeGRangesFromDataFrame(unres, ignore.strand = T, keep.extra.columns = T)
unres.color = "#b20000"


# read in assembly 
file = "~/Desktop/data/genomeWide/Mitchell_CHM1/segdups/CHM1_Falcon_P6_27.fasta.bed5"
asm = read.table(file, header=F)
names(asm) = c(bedform, "conitg", "mapq")
cond = (asm$end - asm$start) > 1000000
asm.long = asm[cond, ]
asm.long = makeGRangesFromDataFrame(asm.long, ignore.strand = T, keep.extra.columns = T)
asm.long.color = "#0000A0"
asm.short = asm[!cond,]
asm.short = makeGRangesFromDataFrame(asm.short, ignore.strand = T, keep.extra.columns = T)
asm.short.color = "#ADD8E6"
asm.short = reduce(asm.short)

# intersect of short asm and unresolved
alreadyDone = GenomicRanges::intersect(unres, asm.short, ignore.strand=T)
alreadyDone = reduce()
alreadyDone.color = "#ADD8E6"


#
# make the plot
#
plotKP <- function(chrs, dopdf=FALSE, name="ideogram", all=FALSE){
  if(dopdf){
    file = paste0("~/Desktop/Public/Mitchell_CHM1/ideogram/", name)
    file = paste0(file, ".pdf")
    pdf(file, width = 30, height = 20)
  }
  if(all){
    kp <- plotKaryotype(genome="hg38", plot.type = 2)
  }else{
    kp <- plotKaryotype(genome="hg38", plot.type = 2, chr = chrs)
  }
  #kpAxis(kp, numticks = 2, tick.pos=c(0.25, 0.75), labels=c("Resolved", "Diverged"), tick.len = 0, side = 2, size = 0.1)
  #kpAxis(kp, numticks = 3, tick.pos=c(0.12,.38,0.80), labels=c("Resolved", "Unresolved", "CHM1 Asm"), tick.len = 0, side = 2, size = 0.1, data.panel = 2)
  kpText(kp, data.panel=1, chr="chr1", x = rep(-15000000,2), y = c(0.25, 0.75), labels = c("Resolved ABP", "Diverged ABP"))
  kpText(kp, data.panel=2, chr="chr1", x = rep(-15000000,2), 
         y = c(0.25,0.75), 
         labels = c("SegDups", "CHM1 Asm"))
  
  
  kpPlotRegions(kp, data=highid, col=highid.color, layer.margin = 0.05, r0=0, r1=0.45)
  kpPlotRegions(kp, data=lowid, col=lowid.color, layer.margin = 0.05,  r0=0.5, r1=1)
  kpPlotRegions(kp, data=alreadyDone, col=alreadyDone.color, layer.margin = 0.05, r0=-0.05, r1=0)
  
  #kpPlotRegions(kp, data=unres, col=unres.color, layer.margin = 0.05, r0=0.0, r1=.45, data.panel = 2)
  kpPlotRegions(kp, data=res, col=res.color, layer.margin = 0.05,  r0=0, r1=0.45, data.panel = 2)
  kpPlotRegions(kp, data=unres, col=unres.color, layer.margin = 0.05, r0=0.0, r1=.45, data.panel = 2)
  
   #kpPlotRegions(kp, data=unres, col=unres.color, layer.margin = 0.02, r0=0.25, r1=.45, data.panel = 2)
  
  kpPlotRegions(kp, data=asm.long, col=asm.long.color, layer.margin = 0.05, r0=.70, r1=0.5, data.panel = 2)
  kpPlotRegions(kp, data=asm.short, col=asm.short.color, layer.margin = 0.05,
               r0=.95, r1=0.75, data.panel = 2, num.layers=2)
  if(dopdf){
    dev.off()
  }
}


for(chr in unique(asm$chr)){
  plotKP(chr, dopdf=T, name=chr)
}

subset = c("chr1","chr2", "chr9", "chr16", "chr14")
plotKP(subset, dopdf=T, name="subset")

plotKP(all=T, dopdf=T)

plotKP(chrs=subset)


