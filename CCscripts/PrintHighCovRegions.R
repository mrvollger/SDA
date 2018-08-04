library(getopt)
options <- matrix(c("cov", "c", 2, "character",
                    "output", "o", 2, "character"), byrow=T, ncol=4)


args <- getopt(options)

#
# Get all coverage
#
cov <- read.table(args$cov)

#
# Splits between contigs.
#
sep <- c(1,which(cov$V1[1:length(cov$V1)-1] != cov$V1[2:length(cov$V1)]), length(cov$V1))

# Names of contigs
names <- unique(cov$V1)



globalSummary <- summary(cov$V4)
globalSdev <- sd(cov$V4)
globalMean <- mean(cov$V4)

GetBlocks <- function(idx){
  pad <- c(-1,idx,max(idx)+2)
  breaks <- which(pad[1:(length(pad)-1)] != pad[2:(length(pad))]-1)
  bi <- which(breaks[2:length(breaks)] - breaks[1:(length(breaks)-1)] > 1)
  starts <- idx[breaks[bi]]
  ends <- idx[breaks[bi+1]-1]
  return(list(s=starts,e=ends))
}


GetHighCovRegions <- function(samp, lim, gap=100) {
  idx <- which(samp >= lim);
  ni <- length(idx)
  if (ni == 0) {
    return(c());
  } else {
    blocks <- GetBlocks(idx)
    nb <- length(blocks$s)
    if (nb <= 1) {
      return(blocks);
    } else {
      
      keep <- rep(T,nb)
      keep[2:nb] = blocks$s[2:nb] - blocks$e[2:nb-1] > gap
      for (i in 2:nb-1) {
        j <- i +1;

        while(j < nb & !keep[j] ) {
          j <- j +1
        }
        blocks$e[i] = blocks$e[j-1]
      }
      blocks <- list(s=blocks$s[keep],e=blocks$e[keep])
      long <- blocks$e - blocks$s > 20      
      return(list(s=blocks$s[long],e=blocks$e[long]))
    }
  }
}

MakeTable <- function(cov, regions, name) {
  if (length(regions$s) == 0) {
    return("");
  }  else {
    return(paste(sapply(1:length(regions$s), function(i) sprintf("%s\t%d\t%d\t%2.2f\t%d\n",name, regions$s[i]*100, regions$e[i]*100, mean(cov[regions$s[i]:regions$e[i]]), regions$e[i]*100- regions$s[i]*100)), sep=""));
  }
}



hcr <- lapply(1:length(names), function(i) GetHighCovRegions(cov$V4[sep[i]:sep[i+1]], globalMean*2- 0.1*globalSdev,1000))


res <- unlist(sapply(1:length(names), function(i) MakeTable(cov$V4[sep[i]:sep[i+1]], hcr[[i]], names[i])))
prunedRes <- res[which(res != "")]
cat(prunedRes, file=args$output,sep="")
