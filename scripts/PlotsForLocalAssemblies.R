#!/usr/bin/env Rscript
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

suppressPackageStartupMessages(library("argparse"))
# create defualt files to run
genome = "Mitchell_CHM13_V2"
genome = "Yoruban"
genome = "Mitchell_CHM1"


res <- Sys.glob(sprintf("~/Desktop/data/genomeWide/%s/segdups/*.mean.resolved", genome) )[1]
unr <- Sys.glob(sprintf("~/Desktop/data/genomeWide/%s/segdups/*.mean.unresolved", genome) )[1]
res <- Sys.glob(sprintf("~/Desktop/data/genomeWide/%s/segdups/*.resolved", genome) )[1]
unr <- Sys.glob(sprintf("~/Desktop/data/genomeWide/%s/segdups/*.unresolved", genome) )[1]
tsv = sprintf("~/Desktop/data/genomeWide/%s/LocalAssemblies/localAssemblyStats.tsv", genome)
des = sprintf("~/Desktop/work/public_html/FigsForABP/%s/", genome)
# create parser object
parser <- ArgumentParser()
parser$add_argument("-t", "--tsv", default=tsv, help="Input tsv file")
parser$add_argument("-o", "--dest", default=des, help="destination folder")
parser$add_argument("-r", "--resolved", default=res, help="list of resolved seg dups")
parser$add_argument("-u", "--unresolved", default=unr, help=" list of unresolved seg dups")
args <- parser$parse_args()
args

#
# set up ggplot theme
#
theme_set(theme_gray(base_size = 18))
h=20
w=30
myTheme =theme(plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
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
names(col4) <- unique(df$Status)
col2=col4[c("Partially Resolved","Resolved")]
col2
myColors = data.frame(Status=names(col4), color = unname(col4))
myColors

#
# function to save plots 
#
mysave <- function(name, p){
  if(args$dest != ""){
    args$dest = paste0(args$dest, "/")
  }
  ggsave(paste0(args$dest, name), plot = p,  width = w, height = h, units = "cm")
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


#
# read in data 
#
#data = "/home/mrvollger/Desktop/data/genomeWide/Mitchell_CHM13_V2/LocalAssemblies/localAssemblyStats.tsv"
#data = "/home/mrvollger/Desktop/data/genomeWide/Mitchell_CHM1/LocalAssemblies/localAssemblyStats.tsv"
#data = "/home/mrvollger/Desktop/data/genomeWide/Yoruban/LocalAssemblies/localAssemblyStats.tsv"
data= args$tsv
df = read.table(data, sep ="\t", header = T)
df$region_in_falcon=as.character(df$region_in_falcon)
df$region_in_falcon <- gsub('\\s+', '', df$region_in_falcon)
df$number_of_CC_groups=as.numeric(df$number_of_CC_groups)
df$copies_in_reference=as.numeric(df$copies_in_reference)
df$numOfPSVs_log10 = log10(df$numOfPSVs)
df$Status=as.character(df$Status)
df[df=="None"]="NA"
df$bestPerID = as.numeric(as.character(df$bestPerID))
df$averagePerID = as.numeric(as.character(df$averagePerID))
df$averageLength = as.numeric(as.character(df$averageLength))
df$bestLength = as.numeric(as.character(df$bestLength))
df[is.na(df)] <- FALSE
df$bestPerID[is.na(df$bestPerID)] <- 0.0

# define the types of resolved 
pr =  "Partially Resolved"
res  = "Resolved"
failed = "Failed"
mAsm = "Multiple Assemblies"
df[df$Status=="Success", "Status" ] = pr 
df[df$Status==pr & df$bestPerID >= 99.8, "Status"] = res
df[df$Status=="MultipleAsm", "Status"] = mAsm

# a poor calcualtions of average reference length 
aveRefLen=c()
for(i in 1:length(df$refRegions)){
  x = as.character(df[i,"refRegions"])
  xx = unlist(strsplit(x,";"))
  xxx = gsub("chr.*:","",xx )
  xxxx = gsub("\\[[0-9]+.[0-9]+\\]","",xxx)
  xxxx 
  lengths = c()
  for( pair in xxxx){
    y = unique(na.omit(as.numeric(unlist(strsplit(pair, "[^0-9]+")))))
    len = y[2] - y[1]
    lengths = c(lengths, len)
  }
  aveRefLen = c(aveRefLen,mean(lengths))
  if(mean(lengths) < 0 ){
    print(aveLen)
    print(x)
    print(xxxx)
  }
}
df$aveRefLen=(aveRefLen)

# add nicer name for num_sam_reads
df$read_count = df$num_of_sam_reads

for( i in 1:length(df$region_in_falcon)){
  region = df[i,]$region_in_falcon
  num = df[i,]$number_of_CC_groups
  num2 = sum(df$region_in_falcon == region)
  if(num!=num2){
    print("hello")
    print(paste0(num,paste0(" ", num2)))
  }
}
df = df[order(df$Status),]
# unique regions in falcon 
dfu = df[!duplicated(df$region_in_falcon),]
# sucesses in the local assemblies 
dfs=df[df$Status=="Resolved" | df$Status=="Partially Resolved", ]
dfpr = dfs[dfs$Status == pr, ]
# sum of each assembly type
sum(df$Status==pr)
sum(df$Status==res)
sum(df$Status==failed)
sum(df$Status==mAsm)

#
# calcualte average percent identity within the reference regions for each collapse 
#
s <- strsplit(as.character(df$refRegions), split = ";")
regions = data.frame(V1 = rep(df$region_in_falcon, sapply(s, length)), V2 = unlist(s))
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
new$region_in_falcon = new$ID
head(new)
dim(df)
df = merge(df, new, by="region_in_falcon")
dim(df)
df = df[order(df$Status),]
df = df[, grep("X[0-9]+", names(df), invert=TRUE, value=T)]
#
# add a column of PSVs per bp in the collapse
#
tail(df)











#
# Counts of each type 
#
p0 <- ggplot(df, aes(Status, fill=Status)) + geom_bar() + 
  labs(y = "Cluster Count") +
  guides(fill=FALSE)+
  scale_fill_manual(values=col4) + myTheme 
p0
mysave("statusHist.png", p0)

#
# number of PSVs per type 
#
p1 <- ggplot(df, aes(Status, numOfPSVs, fill=Status)) + 
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + 
  coord_cartesian(ylim=c(0, 200)) + 
  scale_fill_manual(values=col4) + myTheme
p1
mysave("psv.png", p1)

#
# number of reads per type 
#
p2 <- ggplot(df, aes(Status, read_count, fill=Status)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + 
  coord_cartesian(ylim=c(0, 800)) +
  scale_fill_manual(values=col4) + myTheme
p2
mysave("reads.png", p2)

#
# number of cc groups per type 
#
p3 <- ggplot(df, aes(Status, number_of_CC_groups, fill=Status)) + 
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_fill_manual(values=col4) + myTheme
p3
mysave("numCC.png", p3)

#
# copies in reference vs status 
#
p3 <- ggplot(df, aes(Status, copies_in_reference, fill=Status)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + 
  scale_fill_manual(values=col4) + myTheme
p3
mysave("copies.png", p3)

#
# copies in reference vs status 
#
p3.2 <- ggplot(df, aes(Status, averageRefPerID, fill=Status)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + 
  scale_fill_manual(values=col4) + myTheme + coord_cartesian(ylim=c(90,100))
p3.2
mysave("averageRefPerID.png", p3.2)


#
# average reference length 
#
p4 <- ggplot(df, aes(Status, aveRefLen, fill=Status)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + 
  scale_fill_manual(values=col4) + myTheme
p4
mysave("aveReflen.png", p4)

#
# density of resovled vs partially resolved 
#
p8 = ggplot(dfs, aes(x=length, fill=Status)) + geom_density(alpha=0.8)  + scale_fill_manual(values=col4) +
  scale_x_continuous(labels = comma, breaks = round(seq(min(df$length), max(df$length), by = 25000),1)) +
  theme(axis.text.y = element_blank()) +
  xlab("Assembly length (bp)") + ylab("Density")+ myTheme
p8
mysave("assemblyLength.png",p8)

#
# correlation between cc number and number in reference 
#
corr = round(cor(dfu$copies_in_reference, dfu$number_of_CC_groups,  method = "pearson", use = "complete.obs"), 2)
print(corr)
p5 <- ggplot(dfu, aes(y=copies_in_reference, x=number_of_CC_groups)) + 
  geom_count() + coord_fixed() + geom_abline(intercept = 0, slope = 1, color="red") + 
  ggtitle(paste("R-squared =",as.character(corr))) +
  myTheme
p5
mysave("cc_copies.png", p5)

#
# set up for plotting the number of assemblies that had a certain perID
#
bestPerID = unique(dfs$bestPerID) # gets the unqiue values of perId
ECDF = ecdf(dfs$bestPerID)(unique(dfs$bestPerID)) # gets the ecdf value of each one of those unique items
temp=data.frame(bestPerID, ECDF)
temp$Status = pr
temp[temp$bestPerID >= 99.8, "Status"] = res
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
mysave("cdf.png",p6)

p7 <- p6 +  coord_cartesian(xlim=c(99.5, 100)) +
  geom_segment(aes(x = 99.8, xend = 99.8, y = 0,  yend = miny), color="darkred") +
  geom_segment(aes(x = 0, xend = 99.8, y = miny,  yend = miny), color="darkred") +
  theme(legend.background = element_rect(size=0.5, linetype="solid", color="black"), 
        legend.title=element_blank())
p7
mysave("cdf_zoomed.png", p7)

#
# figure out the number of bases at any given percent identity
#
ID = unique(dfs$bestPerID)
counts = rep(0, length(ID))
names(counts)=ID
for(perID in names(counts)){
  x = df[ df$bestPerID == perID, "length" ] 
  counts[perID] = counts[perID] + sum(x)
}
bestPerID=as.numeric(rep(names(counts), counts))
stopifnot( length(bestPerID) == sum(counts))
BPbyID = data.frame(bestPerID)
BPbyID$Status = pr
BPbyID[BPbyID$bestPerID >= 99.8, "Status"] = res

dfr = dfs[dfs$Status == res,]
sum(dfs[dfs$Status == res, "length"])
sum(dfs$length) 
mean(dfs$length)
median(dfs$length)
mean(dfr$length)
median(dfs[dfs$Status == res, "length"])

bestPerID = unique(BPbyID$bestPerID) # gets the unqiue values of perId
ECDF = ecdf(BPbyID$bestPerID)(unique(BPbyID$bestPerID)) # gets the ecdf value of each one of those unique items
temp=data.frame(bestPerID, ECDF)
temp$Status = pr
temp[temp$bestPerID >= 99.8, "Status"] = res
# adds in an extra row that gits rid of white space between the two fills 
miny=min(temp[temp$Status==res, "ECDF"])
temp2=data.frame(bestPerID=99.8, ECDF=miny , Status=pr)
temp = rbind(temp,temp2)

p6.3 = ggplot() + geom_ribbon(data=temp, aes(x=bestPerID, ymin=0, ymax=ECDF, fill=Status), alpha=0.8 ) + 
  scale_fill_manual(values=col4) + myTheme +
  stat_ecdf(data=BPbyID, aes(bestPerID), pad=FALSE) +
  coord_cartesian(xlim=c(95, 100)) + 
  scale_y_continuous(breaks = seq(0, 1, .1) ) +
  ylab("Fraction of Base Pairs") + xlab("Best Percent Identity Match") +
  theme(legend.background = element_rect(size=0.5, linetype="solid", color="black"), 
        legend.title=element_blank())

#p6.3
#mysave("cdf_by_bp.png", p6.3)

p6.4 <- p6.3 +  coord_cartesian(xlim=c(99.5, 100)) +
  geom_segment(aes(x = 99.8, xend = 99.8, y = 0,  yend = miny), color="darkred") +
  geom_segment(aes(x = 0, xend = 99.8, y = miny,  yend = miny), color="darkred") +
  theme(legend.background = element_rect(size=0.5, linetype="solid", color="black"), 
        legend.title=element_blank())

#p6.4
#mysave("cdf_by_bp_zoomed.png", p6.4)



#
# read in marks data for plot
#
new = df[df$Status=="Resolved", c("aveRefLen", "averagePerID","BestRegionInTheHumanReference")]
new$Status=rep("new", length(new$aveRefLen))
new$ID=seq(1, length(new$aveRefLen))

head(new)
new$averagePerID[new$averagePerID >= 100.0] = 99.999999
dim(new)
x <- new$BestRegionInTheHumanReference

bestRegion = as.data.frame(str_match(x, "^(.*):(.*)-(.*)$")[,-1])
colnames(bestRegion) = c("chr", "start", "end")
bestRegion$chr=(as.character(bestRegion$chr))
bestRegion$start=as.numeric(as.character(bestRegion$start))
bestRegion$end=as.numeric(as.character(bestRegion$end))
new =cbind( new, bestRegion )
head(new)




overlap = perIDfromSegdups(segdups, new)
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
}

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
gw=gw[gw$averagePerID >= 90,]
gw=gw[gw$aveRefLen >= 1000,]

summary(gw)

# add my data to the plot and plot
black <- "#000000"
red <- "#FF0000"
green <- "#228B22"
myColors = c(green, black, red)
#myColors = brewer.pal(3, "Set1")
names(myColors) <- levels(as.factor(gw$Status))

colScale <- scale_colour_manual(name = "Status",values = myColors)
xkb=c("1kb", "10kb", "100kb", "1,000kb")
xbreaks=c(1000, 10000, 100000, 1000000)
ybreaks=seq(90,100,1)
p1 <- ggplot(gw, aes(x=aveRefLen, y=scaled, color=Status) ) + geom_point(size=2) + colScale +
  scale_x_continuous(trans='log10',labels = xkb, breaks = xbreaks) +
  scale_y_continuous(breaks = ybreaks) + 
  xlab("Segmental duplication length (bp)") + ylab("Percent sequence identity") + myTheme 
p1
mysave("newResolvedByPoint.png",p1)


p2 <- ggplot(gw, aes(x=aveRefLen, y=averagePerID, color=Status) ) + geom_density2d(size=3, aes(alpha=..level..)) + guides(alpha=F) +
  colScale +
  scale_x_continuous(trans='log10',labels = xkb, breaks = xbreaks) +
  scale_y_continuous(breaks = ybreaks) + 
  xlab("Segmental duplication length (bp)") + ylab("Percent sequence identity") + myTheme 
p2
mysave("newResolvedByDensity.png", p2)


gw2 = gw[gw$Status != "new",]
p1 <- ggplot(gw2, aes(x=aveRefLen, y=scaled, color=Status) ) + geom_point(size=2) + colScale +
  scale_x_continuous(trans='log10',labels = xkb, breaks = xbreaks) +
  scale_y_continuous(breaks = ybreaks) + 
  xlab("Segmental duplication length (bp)") + ylab("Percent sequence identity") + myTheme 
p1
mysave("unresolvedByPoint.png",p1)

p2 <- ggplot(gw2, aes(x=aveRefLen, y=averagePerID, color=Status) ) + geom_density2d(size=3, aes(alpha=..level..)) + guides(alpha=F) +
  colScale +
  scale_x_continuous(trans='log10',labels = xkb, breaks = xbreaks) +
  scale_y_continuous(breaks = ybreaks) + 
  xlab("Segmental duplication length (bp)") + ylab("Percent sequence identity") + myTheme 
p2
mysave("unresolvedByDensity.png", p2)









#
# quanity failurs per gorup and sucesses per group 
#
region_in_falcon=unique(df$region_in_falcon)
numResolved = c()
numPartiallyResolved = c()
numMultipleAsm = c()
numFailed = c()
for(id in unique(df$region_in_falcon)){
  temp = df[df$region_in_falcon==id,]
  numResolved = c( numResolved, sum(temp$Status == "Resolved")) 
  numPartiallyResolved =c( numPartiallyResolved, sum(temp$Status == "Partially Resolved"))
  numMultipleAsm = c(numMultipleAsm, sum(temp$Status == "Multiple Assemblies"))
  numFailed = c(numFailed, sum(temp$Status == "Failed"))
}
numEach = data.frame(region_in_falcon, numResolved, numPartiallyResolved, numMultipleAsm, numFailed)
head(numEach)

withCounts = merge(dfu, numEach, by="region_in_falcon")
withCounts = withCounts[,c("region_in_falcon",
                           "number_of_CC_groups","copies_in_reference",
                           "numMultipleAsm","numFailed","numPartiallyResolved","numResolved")]
withCounts = merge(withCounts, new, by="region_in_falcon")
head(withCounts)
dim(withCounts)
long <- melt(withCounts, id.vars = c("region_in_falcon", "copies_in_reference", 
                                     "number_of_CC_groups","averagePerID"))
head(long)
p5 <- ggplot(long, aes(y=value, x=copies_in_reference)) + 
  #geom_count()  +  geom_smooth(method='lm',formula=y~x)+
  facet_wrap( ~variable )   +# coord_cartesian(xlim=c(0, 20), ylim=c(0,20)) + 
  myTheme+ stat_summary(fun.y="mean", geom="bar", alpha=0.5) + theme(aspect.ratio = 1) # try with and without
p5
p5 <- ggplot(long, aes(y=value, x=number_of_CC_groups)) + 
  #geom_count()  +  geom_smooth(method='lm',formula=y~x)+
  facet_wrap( ~variable )   + coord_cartesian(xlim=c(0, 20), ylim=c(0,20)) + 
  myTheme+ stat_summary(fun.y="mean", geom="bar", alpha=0.5) + theme(aspect.ratio = 1) # try with and without
p5
p5 <- ggplot(long, aes(y=value, x=averagePerID)) + 
  geom_count() +
  facet_wrap( ~variable )  + coord_cartesian(xlim=c(90, 100)) + 
  myTheme+ stat_summary(fun.y="mean", geom="bar", alpha=0.5) + theme(aspect.ratio = 1) # try with and without
p5
