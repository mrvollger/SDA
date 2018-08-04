#!/usr/bin/env python


import argparse
import ABPUtils
import numpy as np
import sys
import networkx as nx
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser(description="Sort by haplotype")
ap.add_argument("mat", help="snv matrix file")
ap.add_argument("--out", help="Output file", default="/dev/stdout")
ap.add_argument("--graph", help="Write graph to this file.", default=None)
ap.add_argument("--alph", help="Alphabet to use, 3 characters: ref,alt,gap", default='.1n')
ap.add_argument("--cov", help="Average coverage.", type=float, default=60)
ap.add_argument("--score_cutoff", help="Prune connections below this score.",type=int, default=15)
#args = ap.parse_args('assembly.consensus.fragments.snv.mat.categorized')

args = ap.parse_args()

alph     = list(args.alph)
mat      = open(args.mat)
outFile  = open(args.out, 'w')
#freqLine = mat.readline()
#freq     = np.array(freqLine.split())
#print freq
gtList = []
groups = {}
index = 0
groupList = []

coverage = []


(gt, readNames, groupList, groups) = ABPUtils.ReadGenotypeMatrix(mat)
altList = []
refList = []

r=alph[0]
a=alph[1]
g=alph[2]


(altList, refList) = ABPUtils.GetRefAltLists(gt)

coverage = np.array([len(np.where(gt[:,i]!=g)[0]) for i in range(0,gt.shape[1])])

ngt = len(gt)
# Compare distance to members in the group

allGroups = np.array(groups.keys())



allScores = []
nScores = []
scoreIndices = []


for i in range(0,ngt):
    innerMat = []
    innerMis = []
    scores = []
    for j in range(0,ngt):
        if (j == i):
            continue
    
        nMatch = len(np.intersect1d(altList[i],altList[j], assume_unique=True))
        nMis   = len(np.intersect1d(altList[i],refList[j], assume_unique=True))+\
          len(np.intersect1d(refList[i],altList[j], assume_unique=True))

        scores.append([nMatch - nMis, j])

    minAlt = 0
    minRef = 0
    scoreMat = np.array(sorted(scores, reverse=True))
    
    
    if (len(altList[i]) > 0 and len(refList[i]) > 0):
        (start,end) = ABPUtils.GetRange(gt[i], g)
        readCov         = coverage[start]

        bestN = min(readCov, int(args.cov))
        readScores = scoreMat[0:bestN,0]

        minScore = min(readScores)
        scores = np.array([minScore]*int(args.cov))
        scores[0:bestN] = readScores
        outFile.write("\t".join([str(x) for x in scores]) + "\n")

        # record this in a matrix
        scoreIndices.append(i)
        nScores.append(bestN)
        allScores.append(scoreMat[0:bestN,:])

def TripleToHex(x):
    return "#{:02x}{:02x}{:02x}{:02x}".format(int(x[0]*255),int(x[1]*255),int(x[2]*255),int(x[3]*255))

def TripleToTuple(x):
    return "{},{},{}".format(int(x[0]*255),int(x[1]*255),int(x[2]*255))

if (args.graph is not None):
    g = nx.Graph()

    nColors = len(groups.keys())
    groupNames = groups.keys()
    cm = plt.get_cmap("Set1")
    colors = [TripleToHex(cm(int(i*float(float(cm.N)/nColors)))) for i in range(0,len(groupNames))]
    print colors
    groupCM = { groupNames[i]: colors[i] for i in range(0,len(colors))}
    print groupCM
    for i in range(0,ngt):
        g.add_node(i, color=groupCM[groupList[i]])

    idx = 0
    for i in range(0,ngt):
        if (idx >= len(scoreIndices)):
            break
        
        if (scoreIndices[idx] != i):
            continue

        for j in range(0,len(allScores[idx])):
            if (allScores[idx][j][0] < args.score_cutoff):
                break
            g.add_edge(i,allScores[idx][j][1])
        idx+=1
    if (args.graph.find("gml") >= 0):
        nx.write_gml(g, args.graph)
    elif (args.graph.find("gexf") >= 0):
        nx.write_gexf(g, args.graph)
    elif (args.graph.find("graphml") >= 0):
        nx.write_graphml(g, args.graph)
    
    
