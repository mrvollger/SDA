#!/usr/bin/env python


import argparse
import ABPUtils
import numpy as np
import sys
import networkx as nx
import pickle


ap = argparse.ArgumentParser(description="Sort by haplotype")
ap.add_argument("mat", help="snv matrix file")

ap.add_argument("--out", help="Output file", default="/dev/stdout")
ap.add_argument("--graph", help="Write graph to this file.", default=None)
ap.add_argument("--alph", help="Alphabet to use, 3 characters: ref,alt,gap", default='.1n')
ap.add_argument("--simCutoff", help="SNV similarity cutoff.", type=int, default=10)
ap.add_argument("--sharedCutoff", help="Nodes must have this many neighbors in common to retain.", type=int, default=5)
ap.add_argument("--mismatch", help="Penalty for mismatch.",type=float,default=1)
ap.add_argument("--cov", help="Average coverage.", type=float, default=40)
ap.add_argument("--readdm", help="Read distance matrix.", default=None)
ap.add_argument("--writedm", help="Write distance matrix (pickle)", default=None)
ap.add_argument("--readgraph", help="Read graph.", default=None)
ap.add_argument("--maxOvp", help="Maximum # edges per read.", default=None,type=int)


#args = ap.parse_args('assembly.consensus.fragments.snv.mat.categorized')

args = ap.parse_args()

# Set the first step to go to.
step = 0
if (args.readdm):
    step = 1
if (args.readgraph):
    step = 2


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

mat = ABPUtils.ReadGenotypeMatrix(mat)
gt = mat['mat']

altList = []
refList = []

#
# Setup graph, possibly allowing skipping of intermediate steps to
# speed up development
#

if (args.readdm is not None):
    dmFile = open(args.readdm)
    dm = pickle.load(dmFile)

if (step <= 1):
    
    if (args.readdm is None):
        print "Storing distance matrix"
        dm = ABPUtils.StoreDistanceMatrix(gt)



    if (args.writedm is not None):
        dmFile = open(args.writedm, 'wb')
        pickle.dump(dm, dmFile, protocol=pickle.HIGHEST_PROTOCOL)
        
if (args.readgraph is not None):
    graph = ABPUtils.ReadGraph(args.readgraph)
if (step <= 2):
    if (args.readgraph is  None):
        print "Storing graph"
        graph = ABPUtils.DistanceMatrixToGraph(dm, nameLabels=mat['readNames'], colorLabels=mat['groupList'], scoreCutoff=args.simCutoff)
        
    
numPruned = 1
nEdges = graph.number_of_edges()
nRemoved = nEdges
while (nRemoved != 0):
    nEdges = graph.number_of_edges()    
    graph = ABPUtils.PruneGraph(graph, minShared=args.sharedCutoff)
    nRemoved = nEdges - graph.number_of_edges()
    print "Pruned " + str(nRemoved)


if (args.maxOvp is not None):
    graph = ABPUtils.TrimEdges(graph, args.maxOvp)
    
if (args.graph is not None):
    ABPUtils.WriteGraph(graph, args.graph, groupList)
