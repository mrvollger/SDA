library(ggplot2)
library(scales)
library(RColorBrewer)
library(plyr)
library(gridExtra)
library(ggrepel)
library(grid)
library(gtable)



file = "~/Desktop/data/genomeWide/Mitchell_CHM1/LocalAssemblies/betterBlasrMap/careAboutGenes.txt"
importantGenes = read.table(file, header=F)
#importantGenes = c("NPIPA5", "CYP4F", "GOLGA8", "GTF2H2", "DEFA5", "HERC2P3", "CXADPR2", "ZNF705G", "BOLA2")
geneRegex = paste(importantGenes$V1, collapse = '|')
geneRegex

m=1
h=15*m
w=25*m
#
# function to save plots 
#
col2=c("#b20000","#000000")  
names(col2) <- (c(FALSE, TRUE))

mysave <- function(name, p){
  ggsave(name, plot = p,  width = w, height = h, units = "cm")
}

file = "~/Desktop/data/genomeWide/Mitchell_CHM1/LocalAssemblies/betterBlasrMap/clonemapped.pd"
mapped = read.table(file, header=T)

mapped$clonePerID = mapped$perID_by_matches 
#mapped$referencePerID = mapped$X4
mapped$length = mapped$query_end - mapped$query_start
#mapped$frac_in_aln = (mapped$query_end - mapped$query_start) / mapped$query_length
mapped$referenceLength = mapped$bestEnd - mapped$bestStart

mapped$validated = FALSE
mapped[mapped$clonePerID >= 99.8 & mapped$referencePerID <= 99.8 & mapped$frac_in_aln >= 0.9,
       ]$validated = TRUE

summary(mapped[mapped$validated, "length"])
sum(mapped$validated)



bot = 95
top = 100

plots = list()
plots[[1]] =  ggplot(mapped) + geom_point(aes(length, clonePerID, color=validated)) 

plots[[2]] = ggplot(mapped) + geom_point(aes(referenceLength, referencePerID, color=validated))
  #coord_cartesian( ylim = c(bot,top)) 

plots[[3]] = ggplot(mapped) + geom_point(aes(referencePerID, clonePerID, color = validated)) +
  coord_cartesian(xlim=c(bot,top), ylim = c(bot,top)) +
  geom_abline(intercept = 0, slope = 1)
  

validated = mapped[mapped$validated,]
file2 = "~/Desktop/data/genomeWide/Mitchell_CHM1/LocalAssemblies/betterBlasrMap/validated.genes.bed"
genes = read.table(file2, header=F)
genes = genes[, c("V4", "V7")]
colnames(genes) = c("query_name", "Gene")
genes = genes[grep("LOC", genes$Gene, invert=TRUE),]
genes <- ddply(genes, .(query_name), summarize,
               Gene=paste(unique(Gene),collapse=", ") )

valGenes = merge(validated, genes, by="query_name", all.x=T)
valGenes = valGenes[order(valGenes$referencePerID), ]
valGenes$ID = NA
valGenes$GeneShow = NA
counter = 1
maxShownGenes = 15
for(i in 1:length(valGenes$Gene)){
  if( counter <= maxShownGenes & !is.na(valGenes[i,]$Gene) ){
    valGenes[i,]$ID = as.character(counter)
    valGenes[i,]$GeneShow = valGenes[i,]$Gene
    counter = counter + 1
  }
  else if(grepl(geneRegex, valGenes[i,]$Gene)){ # prints if evan thinks it is an important gene
    print(valGenes[i,]$Gene)
    valGenes[i,]$ID = as.character(counter)
    valGenes[i,]$GeneShow = valGenes[i,]$Gene
    counter = counter + 1
  }
}

valGenes$ID
mysize = 3
plots[[4]] = ggplot(valGenes) +
  geom_segment(aes(x=referenceLength, xend=length, y=referencePerID, yend=clonePerID), alpha = 0.25, color="blue") +
  geom_point(aes(x=referenceLength, y=referencePerID), color=col2[[1]], size=mysize) +
  geom_point(aes(x=length, y=clonePerID), color=col2[[2]], size=mysize) +
  geom_label_repel(aes(x=referenceLength, y=referencePerID, label=ID), 
                  size=mysize, color="black", fontface="bold", 
                  arrow = arrow(length = unit(0.005, 'npc')),
                  point.padding = .25,
                  box.padding = .5,
                  #ylim=c(NA, 99.8),
                  segment.size = 0.5) +
  xlab("Alignmnet Length (bp)") +
  ylab("Percent Identity") +
  coord_cartesian(ylim = c(bot+2.3,top)) 


for(i in 1:length(plots)){
  plots[[i]] = plots[[i]] + theme_classic() + 
    scale_x_continuous(labels = comma) + 
    scale_color_manual(values=col2) +
    theme(plot.margin=unit(c(1,1,1,1),"cm"))
  fname = paste0("~/Desktop/Public/Mitchell_CHM1/diverged/", paste0(i, ".all.pdf"))
  mysave(fname, plots[[i]])
}


show = valGenes[!is.na(valGenes$ID), c("ID", "Gene", "clonePerID")]
show$clonePerID = format(round(show$clonePerID, 2), nsmall = 2)
names(show) <- c("Label", "Gene(s)", "% ID")
# Set theme to allow for plotmath expressions

mycolors = c( rep("#FFFFFF", maxShownGenes), rep("#e9e9e9", counter - maxShownGenes))
tsize = 0.3
tt <- gridExtra::ttheme_minimal(
    core = list(fg_params=list(cex = tsize, parse=T), 
                bg_params = list(fill=mycolors  )),
    colhead = list(fg_params=list(cex = tsize, parse=T)))#,
    #rowhead = list(fg_params=list(cex = tsize)))
#tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))


tbl <- tableGrob(show, rows=NULL, theme=tt)
# outline
tbl = gtable_add_grob(tbl, 
                             grobs = rectGrob(gp=gpar(fill=NA, lwd=2)), 
                             t = 1, b = nrow(tbl), l = 1, r = ncol(tbl) )
# header
tbl <- gtable_add_grob(tbl,
                     grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                     t = 1, l = 1, r = ncol(tbl))
# vertical bar
tbl <- gtable_add_grob(tbl,
                       grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                       t = 1, b = nrow(tbl), l = 1)
# highlight
tbl <- gtable_add_grob(tbl,
                       grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                       t = 1, b = nrow(tbl), l = 1)


# Plot chart and table into one object
x = grid.arrange(plots[[4]], tbl,
             ncol=2, 
             widths = c(4,1),
             padding=unit(1, "cm")) 


mysave("~/Desktop/Public/Mitchell_CHM1/diverged/final.pdf", x)

dim(validated)


