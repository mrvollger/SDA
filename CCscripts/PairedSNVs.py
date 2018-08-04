#!/usr/bin/env python


import argparse
import ABPUtils
import numpy as np
import sys
import networkx as nx
import pickle
import intervaltree
#import IPython

ap = argparse.ArgumentParser(description="Sort by haplotype")
ap.add_argument("mat", help="snv matrix file")
ap.add_argument("--minCov", help="Coverage cutoff, to supplement allele fraction", default=10, type=int)
ap.add_argument("--maxCov", help="Coverage cutoff, to supplement allele fraction", default=60, type=int)
ap.add_argument("--vcf", help="Original vcf")
ap.add_argument("--mi",help="Write mutual information to a file", default=None)
ap.add_argument("--graph", help="Write file here", default=None)
ap.add_argument("--adj", help="Write adjacencies here.", default=None)
ap.add_argument("--counts", help="Write ref/alt counts", default=None)
ap.add_argument("--minLRT", help="Minimum log lilelihood ratio test", default=6, type=float)
ap.add_argument("--minNShared",help="Filter min n shared", default=0,type=int)
ap.add_argument("--groups",help="Write groups to this file.", default=None)
ap.add_argument("--cov", help="Estimate of coverage", default=None, type=int)
ap.add_argument("--filt", help="Filt file of snv rates. Used to look up snv position.", default=None)
#args = ap.parse_args('assembly.consensus.fragments.snv.mat.categorized')

args = ap.parse_args()
matFile = open(args.mat)
mat = ABPUtils.ReadGenotypeMatrix(matFile)
gt = mat['mat']

if (args.groups is not None):
    ABPUtils.WriteMinorIndexTable(gt, args.c, mat['groupList'], args.groups)

vcf = []
if (args.vcf is not None):
    vcfFile = open(args.vcf)
    for line in vcfFile:
        if (line[0] != '#'):
            v = line.split()
            vcf.append(v[1])

filtMat = None
if (args.filt is not None):
    filtFile = open(args.filt)
    filtMat = [ [int(i) for i in line.split()[1:]] for line in filtFile]

mi = ABPUtils.FindMutualInformation(gt, args.minCov, args.maxCov, mat['groupList'], miOutFileName=args.mi, minNShared=args.minNShared, filt=filtMat, vcf=vcf, minLRT=args.minLRT)
#
# Determine the frequency of the minor allele
#
(ref,alt) = ABPUtils.GetRefAltCounts(gt)
frac = [0]*len(ref)
for i in range(0,len(ref)):
    denom = ref[i]+alt[i]
    if (denom > 0):
        frac[i] = min(ref[i]/float(denom), alt[i]/float(denom))
    

if (args.counts is not None):
    countsFile = open(args.counts, 'w')
    for i in range(0,len(ref)):
        countsFile.write("{}\t{}\t{}\t{:2.2f}\n".format(i, ref[i], alt[i], frac[i]))
    countsFile.close()
    
nUninformative = 0
if (args.adj is not None):
    adjFile = open(args.adj, 'w')

groupEnum = {}

def NumMinorShared(pair):
    return pair[1]

def NumMajorShared(pair):
    return pair[2]

mig = nx.Graph()

def PassesFilter(mi, args):
    return mi.lrt > args.minLRT

for i in range(0,len(ref)):
    (minor,major) = ABPUtils.GetMinorAllele(ref[i], alt[i], args.minCov, args.maxCov)
    alleleVal=2
    if minor == '.':
        alleleVal = 0
    elif minor == '1':
        alleleVal = 1

    p = 0
    if (len(vcf) > 0):
        p = vcf[i]
    mig.add_node(i,color=0,index=i,allele=alleleVal,pos=p)

groupEnum = {}
for i in range(0,len(ref)):
    if (len(mi[i]) == 0):
        nUninformative+=1
    for j in mi[i].keys():
        if (PassesFilter(mi[i][j], args)):
            e = mi[i][j].maxGroup
            
            if (args.adj is not None):
                if (e not in groupEnum):
                    n = len(groupEnum) +1
                    groupEnum[e] = n
                adjFile.write("{}\t{}\t{}\n".format(i,j,e))
            mig.add_edge(i,j,weight=mi[i][j].totalOverlapping)
            mig[i][j]['color']    = e
            mig[i][j]['weight']   = mi[i][j].fracIMinor

sys.stderr.write("No pairing for " + str(nUninformative) + " sites.\n")

isolatedNodes = []
for n in mig.nodes():
    if len(mig[n].keys()) == 0:
        isolatedNodes.append(n)
mig.remove_nodes_from(isolatedNodes)

if (args.graph is not None):
    ABPUtils.WriteGraph(mig, args.graph)

    
#IPython.embed()
