#!/usr/bin/env python

import argparse
import networkx as nx
import ABPUtils

ap = argparse.ArgumentParser(description="Prune non-solution edges from a graph")
ap.add_argument("--graph", help="Input graph", required=True)
ap.add_argument("--sol", help="Solution", required=True)
ap.add_argument("--out", help="Output file.", required=True)
ap.add_argument("--no-repulsion", help="No repulsion edges", dest='no_repulsion', action='store_true', default=False)
ap.add_argument("--min-value", help="Min value of variable.", default=0.999, type=float, dest="minValue")
args = ap.parse_args()

g = ABPUtils.ReadGraph(args.graph)

#
# Parse cplex solution
#



parseEdges = False
sol = open(args.sol)
subgraph = {}
toRemove = []
nAttract = 0
nRepulse = 0
nRetained = 0
for line in sol:
    if (line.find("Variable Name") >= 0):
        parseEdges = True
        continue
    if (parseEdges == False):
        continue
    var = line.split()[0]
    if (var[0] != "x"):
        continue
    solution = float(line.split()[1])
    edge = var.split("_")
    i = int(edge[1])
    j = int(edge[2])
    if (solution > args.minValue):
        if i not in g or j not in g[i]:
            print "ERROR! Trying to remove edge " + var + " that does not exist."
            continue
            
        toRemove.append((i,j))
        if 'cost' in g[i][j] and g[i][j]['cost'] < 0:
            nRepulse +=1
        else:
            nAttract +=1
        
    else:
        #
        # This pair of edges is fine, record this part of the subgraph
        #
        if (i not in subgraph):
            subgraph[i] = {}
        subgraph[i][j] = True
        nRetained+=1

#
# Now remove everything not listed in the solution
#
print "Removing " + str(nAttract) + " attraction, " + str(nRepulse) + " repulsion from score."
g.remove_edges_from(toRemove)
#print "Should keep " + str(nRetained ) + " in graph."
#toRemove = []
##import pdb
##pdb.set_trace()
#nRemoved = 0
#nKept = 0
#for i in g.nodes():
#    for j in g[i]:
#        if (j > i):
#            if (i not in subgraph or j not in subgraph[i]):
#                if 'cost' in g[i][j] and g[i][j]['cost'] < 0:
#                    nRepulse +=1
#                else:
#                    nAttract +=1
#                toRemove.append((i,j))
#                nRemoved +=1
#            else:
#                nKept +=1
                    
#g.remove_edges_from(toRemove)
#print "Retained " + str(len(g.edges())) + " edges."
import random
toRemove = []
if (args.no_repulsion == True):
    for ep in g.edges():
        e = g[ep[0]][ep[1]]
        if 'cost' in e and e["cost"] < 0:
            toRemove.append(ep)
    print "removing " + str(len(toRemove)) + " repulsion edges"
    g.remove_edges_from(toRemove)
ABPUtils.WriteGraph(g, args.out)


    
