#!/usr/bin/env python

import ABPUtils
import sys
import argparse
import numpy as np

ap = argparse.ArgumentParser(description="Filter the number of neighbors of a source that are shared with neighbors of dest")
ap.add_argument("--graph", help="Input graph")
ap.add_argument("--summary", help="Output summary table.", default=None)
ap.add_argument("--filter", help="Filter edges with this similarity.", type=int, default=0)
ap.add_argument("--filter-graph", help="Write filtered graph here.", dest="filterGraph", default=None)
args = ap.parse_args()

g = ABPUtils.ReadGraph(args.graph)
toRemove = []
if args.summary is not None:
    outFile = open(args.summary,'w')

edgeIndex = 0
nEdges = len(g.edges())
for e in g.edges():
    src = e[0]
    dest = e[1]
    if dest < src:
        continue
    edgeIndex +=1
    if edgeIndex % 100 == 0:
        sys.stderr.write("Processed {}/{} edges\n".format(edgeIndex, nEdges))
    srcNeighbors  = sorted(g[src].keys())
    destNeighbors = sorted(g[dest].keys())
    nShared = len(np.intersect1d(srcNeighbors, destNeighbors, assume_unique=True))
    if args.summary is not None:
        outFile.write("{}\t{}\t{}\t{}\t{}\n".format(src, dest, len(srcNeighbors), len(destNeighbors), nShared))
    if nShared < args.filter:
        toRemove.append(e)

g.remove_edges_from(toRemove)

sys.stderr.write("removing " + str(len(toRemove)) + "\n")
if args.filterGraph is not None:
    ABPUtils.WriteGraph(g, args.filterGraph)
                  
    
