#!/usr/bin/env python
import ABPUtils
import sys
import networkx as nx
g = ABPUtils.ReadGraph(sys.argv[1])
r = []
for e in g.edges():
	if g[e[0]][e[1]]['cost'] < 0:
		r.append(e)

g.remove_edges_from(r)
ABPUtils.WriteGraph(g, sys.argv[2])
 
