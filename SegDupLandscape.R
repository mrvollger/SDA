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
library(data.table)
suppressPackageStartupMessages(library("argparse"))


res <- Sys.glob(sprintf("~/Desktop/data/genomeWide/%s/segdups/*.mean.resolved", genome) )[1]
unr <- Sys.glob(sprintf("~/Desktop/data/genomeWide/%s/segdups/*.mean.unresolved", genome) )[1]
res <- Sys.glob(sprintf("~/Desktop/data/genomeWide/%s/segdups/*.resolved", genome) )[1]
unr <- Sys.glob(sprintf("~/Desktop/data/genomeWide/%s/segdups/*.unresolved", genome) )[1]
euch <- Sys.glob(sprintf("~/Desktop/data/genomeWide/Mitchell_CHM1/LocalAssemblies/euchromatic.hg38.bed", genome) )[1]
tsv = sprintf("~/Desktop/data/genomeWide/%s/LocalAssemblies/localAssemblyStats.tsv", genome)
des = sprintf("~/Desktop/Public/%s/", genome)
# create parser object
parser <- ArgumentParser()
parser$add_argument("-t", "--tsv", default=tsv, help="Input tsv file")
parser$add_argument("-o", "--dest", default=des, help="destination folder")
parser$add_argument("-r", "--resolved", default=res, help="list of resolved seg dups")
parser$add_argument("-u", "--unresolved", default=unr, help=" list of unresolved seg dups")
args <- parser$parse_args()
args
