library(ggplot2)
require(gridExtra)
require(scales)
library(RColorBrewer)

#
# set up ggplot theme
#

# define the types of resolved 
pr =  "Diverged"
res  = "Assembled"
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

myColors


