#!/usr/bin/env python

import ABPUtils
import sys
import networkx as nx
g = ABPUtils.ReadGraph(sys.argv[1])
r = []

for s in g.nodes():
    nAttract = 0
    nRepulse = 0
    for d in sorted(g[s].keys()):
        if g[s][d]['cost'] < 0:
            nRepulse+=1
        else:
            nAttract+=1
    print "{}\t{}\t{}".format(s,nAttract,nRepulse)
 
