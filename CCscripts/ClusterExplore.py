#!/usr/bin/env python

import argparse
import networkx as nx
import ABPUtils
import sys

ap = argparse.ArgumentParser(description="Prune non-solution edges from a graph")
ap.add_argument("--graph", help="Input graph", required=True)
ap.add_argument("--degree-summary", help="Write degree summary", dest="degreeSummary", action='store_true', default=False)
ap.add_argument("--cliques", help="Print clique size for each node.", action='store_true', default=False)
ap.add_argument("--shortest-paths", help="Run shortest path analysis.",
                dest="shortestPaths", default=False, action='store_true')
ap.add_argument("--shortest-paths-file", help="Write shortest paths here.", dest="shortestPathsFile", default=None)

#ap.add_argument("--sol", help="Solution", required=True)
#ap.add_argument("--out", help="Output file.", required=True)

args = ap.parse_args()

g = ABPUtils.ReadGraph(args.graph)

# Now make a version of the graph that does not have any repulsion edges
repl = []
ga = g.copy()
for e in ga.edges():
    if ga[e[0]][e[1]]['cost'] < 0:
        repl.append(e)

ga.remove_edges_from(repl)

if args.cliques:
    for n in ga.nodes():
        print "clique: " + str(n) + "\t" + str(len(ga[n])) + "\t" + str(nx.node_clique_number(ga, n))
    
if args.degreeSummary:
    for s in g:
        nGood = 0
        nBad  = 0
        nNeutral = 0
        p = nx.shortest_path(g, s)
        pa = nx.shortet_path(ga, s)
        for d in g[s]:
            if g[s][d]['cost'] > 0:
                nGood+=1
            if g[s][d]['cost'] == 0:
                nNeutral += 1
            if g[s][d]['cost'] < 0:
                nBad +=1
        if 'color' in g.node[s]:
            col = g.node[s]['color']
        else:
            col = -1
        print str(len(g[s])) + "\t" + str(nGood) + "\t" + str(nBad) + "\t" + str(col)


if args.shortestPaths:
    sys.stderr.write("sp g (all)\n")
    sp = nx.shortest_path(g)
    sys.stderr.write("sp g (no repulsion)\n")    
    spa = nx.shortest_path(ga)
    for s in g:
        for v in sp:
            if v in spa[s] and v in sp[s]:
                import pdb
                pdb.set_trace()
                
                lr = len(sp[s][v])
                la = len(spa[s][v])
                if 'color' in g.node[s] and 'color' in g.node[v]:
                    label = "\t" + str(g.node[s]['color']) + "\t" + str(g.node[v]['color'])
                else:
                    label = ""
                                                                    
                print str(s) + "\t" + str(v) + "\t" + str(lr - la) + label
            
    
