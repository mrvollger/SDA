library(ggplot2)
library(scales)
file = "~/Desktop/data/genomeWide/Mitchell_CHM1/LocalAssemblies/clones/clonemapped.pd"

h=10
w=15
#
# function to save plots 
#
mysave <- function(name, p){
  ggsave(name, plot = p,  width = w, height = h, units = "cm")
}

validated = read.table(file, header=T)
validated$clonePerID = validated$perID_by_matches 
#validated$referencePerID = validated$X4
validated$length = validated$query_end - validated$query_start
validated$validated = FALSE
validated[validated$clonePerID >= 99.8 & validated$referencePerID <= 99.8, ]$validated = TRUE
validated$referenceLength = validated$bestEnd - validated$bestStart

summary(validated[validated$validated, "length"])
validated[which.max(validated$length),]
sum(validated$validated)



bot = 95
top = 100

plots = list()
plots[[1]] =  ggplot(validated) + geom_point(aes(length, clonePerID, color=validated)) 

plots[[2]] = ggplot(validated) + geom_point(aes(referenceLength, referencePerID, color=validated))
  #coord_cartesian( ylim = c(bot,top)) 

plots[[3]] = ggplot(validated) + geom_point(aes(referencePerID, clonePerID, color = validated)) +
  coord_cartesian(xlim=c(bot,top), ylim = c(bot,top)) +
  geom_abline(intercept = 0, slope = 1)
  


for(i in 1:length(plots)){
  plots[[i]] = plots[[i]] + theme_classic() + scale_x_continuous(labels = comma) 
  fname = paste0("~/Desktop/Public/Mitchell_CHM1/diverged/", paste0(i, ".png"))
  mysave(fname, plots[[i]])
}

plots[[1]]
plots[[2]]
plots[[3]]




