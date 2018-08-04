#!/usr/bin/env python

import argparse
import networkx as nx
import ABPUtils

ap = argparse.ArgumentParser(description="Prune non-solution edges from a graph")
ap.add_argument("--graph", help="Input graph", required=True)
ap.add_argument("--sol", help="Solution", required=True)
ap.add_argument("--out", help="Output file.", required=True)
ap.add_argument("--no-repulsion", help="No repulsion edges", dest='no_repulsion', action='store_true', default=False)
ap.add_argument("--min", help="Min value of variable.", default=0.999, type=float)
args = ap.parse_args()

g = ABPUtils.ReadGraph(args.graph)

#
# Parse cplex solution
#
import random

def Vol(g, x, p, b, r):
    w = 0
    for v in b:
        for d in g[v]:
            if d in b:
                w += p[v][d] * x[v][d]
    # always double counting here, cut back to w/2
    return w/2

def Cut(g, S, p):
    c = 0
    for s in S:
        for d in g[s]:
            if d in S:
                w += p[s][d]
    return w/2


def DemaineImmorlicaRound(g, c):
    G = { n: True for n in g.nodes() }

    while len(G) > 0:
        n = g.keys()[0]
        r = 0
        minWeight = 0
        minWeightDest = None
        for d in g[n]:
        
            if w > 0 and w < minWeight:
                minWeight = w
                minWeightDest = d
                
            
        
        
    


parseEdges = False
sol = open(args.sol)
subgraph = {}
toRemove = []
nAttract = 0
nRepulse = 0
for e in g.edges():
    g[e[0]][e[1]]['weight'] = 0

#
# Read the LP to add weights
#

x = {}
p = {}
m = {}

for s in g.nodes():
    x[s] = {}
    m[s] = {}
    p[s] = {}
    for d in g[s]:
        # default to not in the solution
        x[s][d] = 0
        c = g[s][d]['cost']
        if c < 0:
            m[s][d] = abs(c)
            
        else:
            m[s][d] = 0
        if c <= 0:
            p[s][d] = 0
        else:
            p[s][d] = abs(c)

for line in sol:
    if (line.find("Variable Name") >= 0):
        parseEdges = True
        continue
    if (parseEdges == False):
        continue
    var = line.split()[0]
    if (var[0] != "x"):
        continue
    val = float(line.split()[1])

    edge = var.split("_")
    i = int(edge[1])
    j = int(edge[2])
    x[i][j] = val



print "Retaining " + str(nAttract) + " attraction, " + str(nRepulse) + " repulsion"
g.remove_edges_from(toRemove)
toRemove = []
import random
if (args.no_repulsion == True):
    for ep in g.edges():
        e = g[ep[0]][ep[1]]
        if "color" in e and e["color"] == 'red':
            toRemove.append(ep)
    print "removing " + str(len(toRemove)) + " repulsion edges"
    g.remove_edges_from(toRemove)
ABPUtils.WriteGraph(g, args.out)


    
