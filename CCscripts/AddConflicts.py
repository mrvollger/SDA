#!/usr/bin/env python


import argparse
import ABPUtils
import numpy as np
import sys
import networkx as nx
import pickle


ap = argparse.ArgumentParser(description="Convert a weighted graph to a METIS file.")
ap.add_argument("graph", help="Input SNV graph, snvs are connected by an edge if they are in phase")
ap.add_argument("vcf", help="Original vcf file.")
ap.add_argument("sub", help="Subgraph vcf")
ap.add_argument("mat", help="Genotype matrix")
ap.add_argument("--minCov", help="Minimum coverage", type=int,default=20)
ap.add_argument("--maxCov", help="Maximum coverage", type=int, default=80)
ap.add_argument("--minNShared", help="Minimum number of shared sites.", type=int, default=3)
ap.add_argument("--subgraph", help="Write the subgraph here.", default=None)
args = ap.parse_args()

graph = ABPUtils.ReadGraph(args.graph)

vcfFile = open(args.vcf)


vcf = {}
header = []
index = 0

for line in vcfFile:
    if (line[0] == "#"):
        if (line.split()[0] != "#CHROM"):
            header.append(line)
            print header
        continue
    #
    # parse vcf lines
    #001936F 529     snp10   C       T       .       PASS    *       GT:GQ   0/1/1/1:100
    #
    v = line.split()
    # most important are the 2nd & 3rd entries
    contig = v[0]
    vcf[int(v[1])] = [v[2],v[3],v[4],v[5],v[6],v[7],v[8], index]
    index +=1



subVCFFile = open(args.sub)
subset = []
posList = []
for line in subVCFFile:
    if (line[0] == "#"):
        continue
    v = line.split()
    pos = int(v[1])
    index = vcf[pos][-1]
    subset.append(index)
    posList.append(int(v[1]))
varPos = np.array(posList)

matFile = open(args.mat)
mat = ABPUtils.ReadGenotypeMatrix(matFile)
gt = mat['mat']
mi = ABPUtils.FindMutualInformation(gt, args.minCov, args.maxCov, mat['groupList'], minNShared=args.minNShared, conflicts=True, indices=subset, vcf=vcf, pos=varPos)
