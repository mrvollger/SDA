#!/usr/bin/env python


import argparse
import networkx as nx
import ABPUtils
import pdb
ap = argparse.ArgumentParser(description="Add repulsion edges within some radius")
ap.add_argument("--graph", help="Graph file.", required=True)
ap.add_argument("--radius", help="Radius for repulsion", type=int, default=1000)
ap.add_argument("--load", help="At most allow deg(v)*load repulsion edges.",type=int, default=2)
ap.add_argument("--out", help="Output graph.")
args = ap.parse_args()


g = ABPUtils.ReadGraph(args.graph)


sites = [ (int(g.node[n]['pos']), n) for n in g.node]
sites = sorted(sites)

i = 0
nRepulsion = 0
nOK = 0

# store the degrees of nodes before repulsion edges are added

degree = [ len(g[n]) for n in g.node ]
nx.set_edge_attributes(g, 'cost', { e : 1 for e in g.edges() })
component = {}
connectedComponents = [[c for c in comp] for comp in nx.connected_components(g)]
for i in range(0,len(connectedComponents)):
    for n in connectedComponents[i]:
        component[n] = i

for i in range(0,len(sites)-1):
    if degree[i] == 0:
        continue

    post = i+1
    
    while (post < len(sites) and sites[post][0] - sites[i][0] < args.radius):
        post+=1
    nCurRepulsion = 0

    for j in range(i,post):
        if nCurRepulsion/float(degree[i])  > args.load:
            break
        if (j == i):
            continue
        
        else:
            #
            # Make sure repulsion edges are not added to disconnected nodes
            #
            if (g.has_edge(sites[j][1], sites[i][1]) == False):
                if (component[sites[i][1]] == component[sites[j][1]]):
                    nCurRepulsion += 1
                    g.add_edge(sites[i][1], sites[j][1], cost=-1, weight=0.1)
            else:
                nOK +=1
    for j in range(post, len(sites)):
        if (g.has_edge(sites[i][1], sites[j][1])):
            nOK +=1
    nRepulsion += nCurRepulsion

print str(nRepulsion) + "\t" + str(nOK)
ABPUtils.WriteGraph(g, args.out)
