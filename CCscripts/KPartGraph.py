#!/usr/bin/env python

import argparse
import ABPUtils
import numpy as np
import sys
import networkx as nx
import matplotlib.pyplot as plt
import pickle
import math

ap = argparse.ArgumentParser(description="Convert a weighted graph to a METIS file.")
ap.add_argument("graph", help="Input graph file.")
ap.add_argument("--shuffle", help="Shuffle initial partitions", default=False, action='store_true')
ap.add_argument("--k", help="Number of partitions", default=None, type=int)
ap.add_argument("--f", help="Fraction of partitions k is total args supplied", nargs="+", type=float, default=None)
ap.add_argument("--out", help="Output file.",default="/dev/stdout")


args = ap.parse_args()

g = ABPUtils.ReadGraph(args.graph)
ABPUtils.RemoveIsolated(g)


edgesWithData = g.edges(data=True)

weights = [e[2]['capacity'] for e in edgesWithData]

edges =g.edges()

nNodes=g.number_of_nodes()

if (args.k is not None):
    partFrac = [1.0/args.k]*args.k
elif (args.f is not None):
    partFrac = args.f

partSizes = np.asarray(partFrac)*nNodes
partSizes = partSizes.astype(np.int32)



partitions = {}
cur = 0

if (len(partSizes) > 1):
    pre = sum(partSizes[0:-1])
    partSizes[-1] = nNodes - pre
else:
    partSizes[0] = nNodes

nodes = g.nodes()
for i in range(0,len(partFrac)):
    for j in range(cur,cur+partSizes[i]):
        partitions[nodes[j]] = i
    cur+=partSizes[i]

import random
if (args.shuffle):
    keys = partitions.keys()
    random.shuffle(keys)
    for i in range(0,len(keys)-1,2):
        t = partitions[keys[i]]
        partitions[keys[i]] = partitions[keys[i+1]]
        partitions[keys[i+1]] = t

vg = ABPUtils.VectorizeGraph(g)
        
ABPUtils.OptimizePartitions(partitions, vg, g.nodes())

optScore = ABPUtils.EvaluateTruthPartitions(g, vg)
print "optScore: " + str(optScore)

nx.set_node_attributes(g, 'part', { i: partitions[i] for i in g.nodes() } )

ABPUtils.WriteGraph(g, args.out)
