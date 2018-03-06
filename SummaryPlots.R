#!/usr/bin/env Rscript
#lib="/home/mrvollger/anaconda3/lib/R/library"
#.libPaths(lib)
#.libPaths( c( .libPaths(), lib ))
library(Cairo) 
library(ggplot2)
library(plyr)
require(gridExtra)
require(scales)
library(RColorBrewer)
library(reshape2)
#install.packages("tidyr")
library(tidyr)
library(dplyr)
library(splitstackshape)
library(stringr)
library(data.table)
#library(networkD3)
library(bedr)
library(evaluate)
suppressPackageStartupMessages(library("argparse"))

# create defualt files to run
genome = "Mitchell_CHM13_V2"
genome = "Yoruban"
genome = "Mitchell_CHM1"
# this is the one with jason chins assembly
genome = "Mitchell_CHM1_V2"


faifile = sprintf("~/Desktop/data/genomeWide/%s/LocalAssemblies/all.ref.fasta.fai", genome)
res <- Sys.glob(sprintf("~/Desktop/data/genomeWide/%s/segdups/*.mean.resolved", genome) )[1]
unr <- Sys.glob(sprintf("~/Desktop/data/genomeWide/%s/segdups/*.mean.unresolved", genome) )[1]
res <- Sys.glob(sprintf("~/Desktop/work/assemblies/CHM1/GCA_001297185.1_PacBioCHM1_r2_GenBank_08312015/Segdups/asm.resolved", genome) )[1]
unr <- Sys.glob(sprintf("~/Desktop/work/assemblies/CHM1/GCA_001297185.1_PacBioCHM1_r2_GenBank_08312015/Segdups/asm.unresolved", genome) )[1]
#euch <- Sys.glob(sprintf("~/Desktop/data/genomeWide/Mitchell_CHM1/LocalAssemblies/euchromatic.hg38.bed", genome) )[1]
tsv = sprintf("~/Desktop/data/genomeWide/%s/LocalAssemblies/localAssemblyStats.tsv", genome)
des = sprintf("~/Desktop/data/genomeWide/%s/plots/", genome)
# create parser object
parser <- ArgumentParser()
parser$add_argument("-t", "--tsv", default=tsv, help="Input tsv file")
parser$add_argument("-o", "--dest", default=des, help="destination folder")
parser$add_argument("-r", "--resolved", default=res, help="list of resolved seg dups")
parser$add_argument("-u", "--unresolved", default=unr, help=" list of unresolved seg dups")
args <- parser$parse_args()
args

if( ! dir.exists(args$dest)){
  dir.create(args$dest)
}


#
# set up ggplot theme
#

# define the types of resolved 
pr =  "Diverged"
res  = "Resolved"
failed = "Failed"
mAsm = "Multiple Assemblies"
theme_set(theme_gray(base_size = 24))
h=20
w=30
font="serif"
font="Comic Sans MS"
font="Times"
myTheme =theme(plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5, family=font),
               text = element_text(),
               panel.background = element_blank(),
               panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),
               plot.background = element_rect(colour = NA),
               axis.title = element_text(face = "bold",size = rel(1)),
               axis.title.y = element_text(angle=90,vjust =2),
               axis.title.x = element_text(vjust = -0.2),
               axis.text = element_text(), 
               axis.line = element_line(colour="black"),
               axis.ticks = element_line(),
               #panel.grid.major = element_line(colour="#f0f0f0"),
               #panel.grid.minor = element_blank(),
               plot.margin=unit(c(10,5,5,5),"mm"),
               #strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
               strip.text = element_text(face="bold")
)
col4=c("#b20000","#660000","#A9A9A9","#000000")  
names(col4) <- (c(failed,mAsm,pr, res))
col4
col2=col4[c(pr, res)]
col2
myColors = data.frame(Status=names(col4), color = unname(col4))
myColors
#
# function to save plots 
#
mysave <- function(name, p){
  if(args$dest != ""){
    args$dest = paste0(args$dest, "/")
    file = paste0(args$dest, name)
    ggsave(file, plot = p,  width = w, height = h, units = "cm")
    # save as an svg as well
    svg(paste0(file, ".svg"), width = w/2, height = h/2)
    print(p)
    dev.off()
  }
}

#
# this funciton calcualtes overlap between two data frames with columns chr, start, end
#
perIDfromSegdups <- function(rangesA, rangesB){
  ranges <- merge(rangesA,rangesB,by="chr",suffixes=c("A","B"))
  for( x in c("startA", "startB", "endA", "endB")){
    ranges[[x]] = as.numeric(as.character(ranges[[x]]))  
  }
  postions =  (ranges$startB > ranges$startA & ranges$startB < ranges$endA ) |
    (ranges$endB > ranges$startA & ranges$endB < ranges$endA ) |
    (ranges$startB < ranges$startA & ranges$endB > ranges$endA ) 
  ranges = ranges[postions,]
  ranges <- transform(ranges, start = pmax(startA, startB))
  ranges <- transform(ranges, end = pmin(endA, endB))
  ranges <- transform(ranges, overlap = end - start)
  return(ranges)
}

#
# this is a function to read in collapsed and uncollapsed segdups
#
readIn <- function(filename, name){
  df=read.table(filename, sep ="\t", header = F)
  colnames(df) = c("chr", "start", "end", "averagePerID")
  df$start=as.numeric(df$start)
  df$end=as.numeric(df$end)
  df$aveRefLen=df$end-df$start
  # must adjust to be between zero and 100
  if( max(df$averagePerID <= 1.0) ){
    df$averagePerID = df$averagePerID * 100
  }
  df$averagePerID[df$averagePerID >= 100.0] = 99.999999
  df$Status=rep(name, length(df$aveRefLen))
  return(df)
}

resolved = readIn(args$resolved,"resolved")
head(resolved)
unresolved = readIn(args$unresolved, "unresolved")
head(unresolved)
segdups = rbind(resolved, unresolved)
head(segdups)

euchRes = readIn(args$resolved,"resolved")
head(resolved)
unresolved = readIn(args$unresolved, "unresolved")
head(unresolved)
segdups = rbind(resolved, unresolved)
head(segdups)


#
# read in data 
#

# remove multiple instances of the header 
lines <- readLines(args$tsv)
header = lines[1]
these = !grepl("reference_name", lines)
lines = lines[ these ] 
lines = c(header, lines)
length(lines)
all = read.table(text=lines, sep ="\t", header = T, stringsAsFactors=FALSE)
setDT(all)
# convert to percentage vs fraciton
all$bestPerID = all$perID_by_matches
sum(all$Status == "Diverged")
sum(all$Status == "Resolved")

#
# calcualte average percent identity within the reference regions for each collapse 
#
s <- strsplit(as.character(all$refRegions), split = ";")
regions = data.frame(V1 = rep(all$collapse, sapply(s, length)), V2 = unlist(s))
head(regions)
split = as.data.frame(str_match(regions$V2, "^(.*):(.*)-([0-9]*)(.*)$")[,-1])
regions = cbind(regions$V1, split)
colnames(regions) = c("ID", "chr", "start", "end", "useless")
regions = regions[,c("ID", "chr", "start", "end")]
head(regions)

overlap = perIDfromSegdups(segdups, regions)
overlap = overlap[!duplicated(overlap), ]
head(overlap)

new = data.frame(ID = unique(regions$ID) )
#new$averageRefLen = rep(0, length(new$ID))
new$averageRefPerID = rep(0, length(new$ID))
for(id in unique(overlap$ID) ){
  temp = overlap[overlap$ID == id,]
  total = sum(temp$overlap)
  #print(temp$averagePerID)
  #print(temp$overlap)
  product = sum( temp$averagePerID * temp$overlap )
  perID = product/total
  #new[new$ID == id, ]$averageRefLen = total
  new[new$ID == id, ]$averageRefPerID = perID
}
new$collapse = new$ID
head(new)
dim(all)
all = merge(all, new, by="collapse")
dim(all)
all = all[order(all$Status),]

#
# add a column of PSVs per bp in the collapse
#
all$PSVsPer1K = all$totalPSVs / all$collapseLen * 1000
all$PSVsPerClusterPer1K = all$numPSVs / all$collapseLen * 1000

all$Status[all$Status == "Partially Resolved"] = pr

# remove multiple assemblies that appear more than once. I want them overall for mapping them, but i do not want them for this summary
dim(all)
df =all[ !duplicated(all, by=c("collapse", "Status", "CC_ID"), fromLast = F) ]
dim(df)
#
# reorder Status factors for ggplot
#
df$Status <- factor(df$Status, levels=c(failed, mAsm, pr, res))
#
# remove some NAs
#
df$Length[is.na(df$Length)] = 0
df$bestPerID[is.na(df$bestPerID)] = 0
# show data tpye of each col
sapply(df, class)
# unique regions in falcon 
dfu = df[!duplicated(df$collapse),]
# sucesses in the local assemblies 
dfs=df[df$Status==res | df$Status==pr, ]
dfpr = dfs[dfs$Status == pr, ]
# sum of each assembly type
counts = c(sum(df$Status==pr), sum(df$Status==res), sum(df$Status==failed), sum(df$Status==mAsm))
names(counts) = c(pr,res,failed,mAsm)
counts
statusCounts = data.frame(Status = c(pr, res, failed, mAsm), 
                     counts = counts)
dim(df)
statusCounts

p0 <- ggplot(df, aes(Status, fill=Status)) + geom_bar() + 
  labs(y = "Cluster Count") +
  guides(fill=FALSE)+
  geom_text(data=statusCounts, aes(y=counts, x=Status, label=counts),vjust=-.1, size = 8) +
  scale_fill_manual(values=col4) + myTheme 
p0
mysave("statusHist.pdf", p0)



nodes <- data.frame(name=c(
"Clusters",
failed,
"Assembled",
mAsm,
"Successful Assemblies",
pr,
res))
nodes

source=c()
target=c()
value=c()
# failed
    source = c(source, 0)
    target = c(target, 1)
    value=c(value, counts[failed])
# mAsm
    source = c(source, 0)
    target = c(target, 2)
    value=c(value, counts[res] + counts[pr] + counts[mAsm] )
    
    source = c(source, 2)
    target = c(target, 3)
    value=c(value, counts[mAsm])
    
# res
    source = c(source, 2)
    target = c(target, 4)
    value=c(value, counts[res] + counts[pr] )
    
    source = c(source, 4)
    target = c(target, 6)
    value=c(value, counts[res])
    
# diverged
    source = c(source, 4)
    target = c(target, 5)
    value=c(value, counts[pr])
    
    
links = data.frame(source, target, value)



#
# density of reference assembly length 
#
fai = read.table(faifile)
names(fai) = c("Collapse","Collapse Length","b","a")
fai[["Collapse Length"]][ fai[["Collapse Length"]] > 150000 ] = 150000
p8 = ggplot(fai, aes(x=fai[["Collapse Length"]] )) + geom_histogram(alpha=0.8, fill="black", binwidth=3000)  + scale_fill_manual(values=col4) +
  scale_x_continuous(labels = comma, breaks = seq(0, 150000, by = 25000) ) +
  xlab("Collapse length (bp)") + ylab("Counts")+ myTheme
p8
mysave("CollapseLength.pdf", p8)


#
# number of PSVs per type 
#
p1 <- ggplot(df, aes(Status, PSVsPer1K, fill=Status)) + 
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + 
  coord_cartesian(ylim=c(0, 100)) + 
  scale_fill_manual(values=col4) + myTheme
p1
mysave("PSVsPer1K.pdf", p1)


#
# number of PSVs per type per cluster 
#
p1 <- ggplot(df, aes(Status, PSVsPerClusterPer1K, fill=Status)) + 
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + 
  coord_cartesian(ylim=c(0, 7)) + 
  scale_fill_manual(values=col4) + myTheme
p1
mysave("PSVsPerClusterPer1K.pdf", p1)
#
# number of PSVs per type 
#
p1 <- ggplot(df, aes(Status, numPSVs, fill=Status)) + 
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + 
  coord_cartesian(ylim=c(0, 100)) + 
  scale_fill_manual(values=col4) + myTheme
p1
mysave("numPSVs.pdf", p1)



#
# number of reads per type 
#
p2 <- ggplot(df, aes(Status, numReads, fill=Status)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + 
  coord_cartesian(ylim=c(0, 800)) +
  scale_fill_manual(values=col4) + myTheme
p2
mysave("reads.pdf", p2)


#
# number of cc groups per type 
#
p3 <- ggplot(df, aes(Status, numOfCCgroups, fill=Status)) + 
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_fill_manual(values=col4) + myTheme
p3
mysave("numOfCCgroups.pdf", p3)

#
# copies in reference vs status 
#
p3 <- ggplot(df, aes(Status, copiesInRef, fill=Status)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + 
  scale_fill_manual(values=col4) + myTheme
p3
mysave("copies.pdf", p3)

#
# copies in reference vs status 
#
p3.2 <- ggplot(df, aes(Status, averageRefPerID, fill=Status)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + 
  scale_fill_manual(values=col4) + myTheme + coord_cartesian(ylim=c(90,100))
p3.2
mysave("averageRefPerID.pdf", p3.2)


#
# average reference length 
#
p4 <- ggplot(df, aes(Status, aveRefLength, fill=Status)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + 
  scale_fill_manual(values=col4) + myTheme
p4
mysave("aveRefLength.pdf", p4)

#
# density of resovled vs partially resolved 
#
p8 = ggplot(dfs, aes(x=Length, fill=Status)) + geom_density(alpha=0.8)  + scale_fill_manual(values=col4) +
  scale_x_continuous(labels = comma, breaks = round(seq(min(df$Length), max(df$Length), by = 25000),1)) +
  theme(axis.text.y = element_blank()) +
  xlab("Assembly length (bp)") + ylab("Density")+ myTheme
p8
mysave("assemblyLength.pdf",p8)

#
# correlation between cc number and number in reference 
#
corr = round(cor(dfu$copiesInRef, dfu$numOfCCgroups), 2)
print(corr)
p5 <- ggplot(dfu, aes(y=copiesInRef, x=numOfCCgroups)) + 
  geom_count() + coord_fixed() + geom_abline(intercept = 0, slope = 1, color="red") + 
  ggtitle(paste("R-squared =",as.character(corr))) +
  myTheme
p5
mysave("cc_copies.pdf", p5)



#
# correlation numer of failuses and number of copies 
#
corr = round(cor(dfu$copiesInRef, dfu$numF), 2)
print(corr)
p5 <- ggplot(dfu, aes(y=copiesInRef, x=numF)) + 
  geom_count() +
  ggtitle(paste("R-squared =",as.character(corr))) +
  myTheme
p5
mysave("failsVsCopies.pdf", p5)
#
# correlation numer of failuses and number of copies 
#
corr = round(cor(dfu$copiesInRef, dfu$numOfCCgroups), 2)
print(corr)
p5 <- ggplot(dfu, aes(y=numOfCCgroups, x=numF)) + 
  geom_count() +
  ggtitle(paste("R-squared =",as.character(corr))) +
  myTheme
p5
mysave("failsVsClusters.pdf", p5)
#
# correlation numer of failuses and number of copies 
#
dfuu = dfu[ dfu$averageRefPerID > 80, ]
corr = round(cor(dfuu$averageRefPerID, dfuu$numOfCCgroups), 2)
print(corr)
p5 <- ggplot(dfu, aes(y=averageRefPerID, x=numF)) + 
  geom_count() +
  ggtitle(paste("R-squared =",as.character(corr))) +
  coord_cartesian(ylim=c(80, 100)) +
  myTheme
p5
mysave("failsVsPerID.pdf", p5)






#
# set up for plotting the number of assemblies that had a certain perID
#
bestPerID = unique(dfs$bestPerID) # gets the unqiue values of perId
ECDF = ecdf(dfs$bestPerID)(unique(dfs$bestPerID)) # gets the ecdf value of each one of those unique items
temp = data.frame(bestPerID, ECDF)
temp$Status = ifelse(temp$bestPerID >= 99.8, res, pr)
# adds in an extra row that gits rid of white space between the two fills 
miny=min(temp[temp$Status==res, "ECDF"])
temp2=data.frame(bestPerID=99.8, ECDF=miny , Status=pr)
temp = rbind(temp,temp2)

p6 = ggplot() + geom_ribbon(data=temp, aes(x=bestPerID, ymin=0, ymax=ECDF, fill=Status), alpha=0.8 ) + 
  scale_fill_manual(values=col4) + myTheme +
  stat_ecdf(data=dfs, aes(bestPerID), pad=FALSE) +
  coord_cartesian(xlim=c(95, 100)) + 
  scale_y_continuous(breaks = seq(0, 1, .1) ) +
  ylab("Fraction of Assemblies") + xlab("Best Percent Identity Match") +
  theme(legend.background = element_rect(size=0.5, linetype="solid", color="black"), 
        legend.title=element_blank())
p6
mysave("cdf.pdf",p6)

p7 <- p6 +  coord_cartesian(xlim=c(99.5, 100)) +
  geom_segment(aes(x = 99.8, xend = 99.8, y = 0,  yend = miny), color="darkred") +
  geom_segment(aes(x = 0, xend = 99.8, y = miny,  yend = miny), color="darkred") +
  theme(legend.background = element_rect(size=0.5, linetype="solid", color="black"), 
        legend.title=element_blank())
p7
mysave("cdf_zoomed.pdf", p7)

#
# figure out the number of bases at any given percent identity
#
ID = unique(dfs$bestPerID)
counts = rep(0, length(ID))
names(counts)=ID
for(perID in names(counts)){
  x = df[ df$bestPerID == perID, "Length" ] 
  counts[perID] = counts[perID] + sum(x)
}
bestPerID=as.numeric(rep(names(counts), counts))
stopifnot( length(bestPerID) == sum(counts))
BPbyID = data.frame(bestPerID)
BPbyID$Status = pr
BPbyID[BPbyID$bestPerID >= 99.8, "Status"] = res

dfr = dfs[dfs$Status == res,]
sum(dfs[dfs$Status == res, "Length"])
sum(dfs$Length) 
mean(dfs$Length)
median(dfs$Length)
mean(dfr$Length)
median(dfs[dfs$Status == res]$Length)
summary(dfr$Length)


bestPerID = unique(BPbyID$bestPerID) # gets the unqiue values of perId
ECDF = ecdf(BPbyID$bestPerID)(unique(BPbyID$bestPerID)) # gets the ecdf value of each one of those unique items
temp=data.frame(bestPerID, ECDF)
temp$Status = pr
temp[temp$bestPerID >= 99.8, "Status"] = res
# adds in an extra row that gits rid of white space between the two fills 
miny=min(temp[temp$Status==res, "ECDF"])
temp2=data.frame(bestPerID=99.8, ECDF=miny , Status=pr)
temp = rbind(temp,temp2)

ymax = round(sum(dfs$Length)/1000000) # megabases of sequences
ymax
p6.3 = ggplot() + geom_ribbon(data=temp, aes(x=bestPerID, ymin=0, ymax=ECDF, fill=Status), alpha=0.8 ) + 
  scale_fill_manual(values=col4) + myTheme +
  stat_ecdf(data=BPbyID, aes(bestPerID), pad=FALSE) +
  coord_cartesian(xlim=c(95, 100)) + 
  scale_y_continuous(breaks = seq(0, 1, 0.1), labels = round(seq(0, ymax, ymax/10)) ) +
  ylab("Mb of Assembly") + xlab("Best Percent Identity Match") +
  theme(legend.background = element_rect(size=0.5, linetype="solid", color="black"), 
        legend.title=element_blank())

p6.4 <- p6.3 +  coord_cartesian(xlim=c(97.5, 100)) +
  geom_segment(aes(x = 99.8, xend = 99.8, y = 0,  yend = miny), color="darkred") +
  geom_segment(aes(x = 0, xend = 99.8, y = miny,  yend = miny), color="darkred") +
  theme(legend.background = element_rect(size=0.5, linetype="solid", color="black"), 
        legend.title=element_blank())
#
# this takes forever and lots of memory, so only plot if I set this if statemtn
#
if(F){
  p6.3
  mysave("cdf_by_bp.pdf", p6.3)
  p6.4
  mysave("cdf_by_bp_zoomed.pdf", p6.4)
}




#
# read in marks data for plot
#
new = df[df$Status=="Resolved", c("bestMatch")]
new$averagePerID = df[df$Status=="Resolved", "averageRefPerID"]
new$aveRefLen = df[df$Status=="Resolved", "aveRefLength"]
new$chr = df[df$Status=="Resolved", "bestChr"]
new$start= df[df$Status=="Resolved", "bestStart"]
new$end = df[df$Status=="Resolved","bestEnd"]
new$Status=rep("ABP", length(new$aveRefLen))
new$ID=seq(1, length(new$aveRefLen))
new$averagePerID[new$averagePerID >= 100.0] = 99.999999
dim(new)
head(new)

#x <- new$BestRegionInTheHumanReference
#bestRegion = as.data.frame(str_match(x, "^(.*):(.*)-(.*)$")[,-1])
#colnames(bestRegion) = c("chr", "start", "end")
#bestRegion$chr=(as.character(bestRegion$chr))
#bestRegion$start=as.numeric(as.character(bestRegion$start))
#bestRegion$end=as.numeric(as.character(bestRegion$end))
#new =cbind( new, bestRegion )
#head(new)


names(segdups)
head(segdups)
overlap = perIDfromSegdups(segdups[, -6], new)
head(overlap)
# remove thigns that have overlap with zero segdups
new = new[new$ID %in% unique(overlap$ID), ]
for(id in unique(overlap$ID) ){
  temp = overlap[overlap$ID == id,]
  total = sum(temp$overlap)
  totalOverlap = sum(temp$overlap)
  if(total > totalOverlap ){
    print("yo overlapping things, lame")
  }
  product = sum( temp$averagePerIDA * temp$overlap )
  perID = product/total
  new[new$ID == id, ]$aveRefLen = totalOverlap
  new[new$ID == id, ]$averagePerID = perID
  print(dim(temp))
  print(id)
}
summary(new)

new$dec=new$averagePerID - floor(new$averagePerID)
new$dec = rescale( new$dec , to=c(2/3,1), from=c(0,1) )
new$scaled = floor(new$averagePerID) + new$dec

resolved$dec=resolved$averagePerID - floor(resolved$averagePerID)
resolved$dec = rescale( resolved$dec , to=c(0,1/3), from=c(0,1) )
resolved$scaled = floor(resolved$averagePerID) + resolved$dec

unresolved$dec=unresolved$averagePerID - floor(unresolved$averagePerID)
unresolved$dec = rescale( unresolved$dec , to=c(1/3,2/3), from=c(0,1) )
unresolved$scaled = floor(unresolved$averagePerID) + unresolved$dec

temp = merge(resolved, unresolved,all=TRUE)
gw = merge(temp, new,all=TRUE)
dim(gw)
gw=gw[gw$averagePerID >= 90,]
dim(gw)
gw=gw[gw$aveRefLen >= 1000,]
dim(gw)
summary(gw)


#
# function for plotting this data
#
plotSegDups<-function(data, name1, name2, myColors){
  #data = gw
  colScale <- scale_colour_manual(name = "Status", values = myColors)
  xkb=c("1", "10", "100", "1,000")
  xbreaks=c(1000, 10000, 100000, 1000000)
  print(length(xbreaks))
  print(length(xkb))
  ybreaks=seq(90,100,1)
  print("zero done")
  p1 <- ggplot(data, aes(x=aveRefLen, y=scaled, color=Status) ) + geom_point(size=2) + colScale +
    scale_x_continuous(trans='log10',labels = xkb, breaks = xbreaks) +
    scale_y_continuous(breaks = ybreaks, labels=ybreaks) + 
    xlab("Segmental duplication length (kb)") + ylab("Percent sequence identity") + myTheme 
  p1
  mysave(name1, p1)
  
  print("one done")
  
  p2 <- ggplot(data, aes(x=aveRefLen, y=averagePerID, color=Status) ) +
    geom_density2d(size=3, aes(alpha=..level..)) + guides(alpha=F) +
    colScale +
    scale_x_continuous(trans='log10',labels = xkb, breaks = xbreaks) +
    scale_y_continuous(breaks = ybreaks) + 
    xlab("Segmental duplication length (kb)") + ylab("Percent sequence identity") + myTheme 
  p2
  mysave(name2, p2)
  
  list(p1,p2)
}

# add my data to the plot and plot
black <- "#000000"
red <- "#FF0000"
green <- "#00b2b2"
myColors = c(green, black, red)
names(myColors) <- levels(as.factor(gw$Status))

plots = ( plotSegDups(gw,"newResolvedByPoint.pdf" ,"newResolvedByDensity.pdf", myColors) )
plots[1]
plots[2]

gw2 = gw[gw$Status != "ABP",]
plots2 = ( plotSegDups(gw2,"ResolvedByPoint.pdf" ,"ResolvedByDensity.pdf", myColors) )
plots2[1]
plots2[2]



library(bedr)
bed = gw[c("chr","start","end", "Status")]
bed$chr = as.character.factor(bed$chr)
is.valid.region(bed)
merged = bedr.merge.region(bed, check.chr = FALSE)
merged$len = merged$end - merged$start
sum(merged$len)

numres = bed[bed$Status != "unresolved", ]
numres = bedr.merge.region(numres, check.chr = FALSE)
numres$len = numres$end - numres$start
sum(numres$len)

unres = bed[bed$Status == "unresolved", ]
unres = bedr.merge.region(unres, check.chr = FALSE)
unres$len = unres$end - unres$start
sum(unres$len)

sum(numres$len) / sum(merged$len)

1-sum(unres$len ) / sum(merged$len) 

#
# Run a random forest to detemine what partially resolved ones might actually be correct 
# 
library(randomForest)

#dff <- as.data.frame(unclass(all))
dff <- as.data.frame(unclass(dfs))
dff$Status = droplevels(dff$Status)
sapply(dff, class)

# model fit
# note that you must turn the ordinal variables into factor or R wont use
# them properly
model <- randomForest(formula = Status ~ totalPSVs + numPSVs + aveRefLength + averageRefPerID + numReads + 
                        Length + copiesInRef + PSVsPerClusterPer1K + collapseLen + numOfCCgroups + totalReads, 
                      data  = dff,
                      ntree=2000, importance=TRUE)
#plot of model accuracy by class
plot(model)
model

# Importance of each predictor.
print(importance(model))
varImpPlot(model)

dff$prediction = predict(model)
candidates = dff[ dff$Status == pr & dff$prediction == res,]
candidates = candidates[order(candidates$bestPerID),]
candidates = candidates[candidates$bestPerID > 96.0 & candidates$bestPerID < 99.0,]
candidates


