#!/usr/bin/env python

import argparse
import networkx
import ABPUtils
ap = argparse.ArgumentParser(description="Print components of file.")
ap.add_argument("operation", help="Operation to run")
ap.add_argument("graph", help="Input graph file.")
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
args = ap.parse_args()

graph = ABPUtils.ReadGraph(args.graph)
outFile = ABPUtils.OpenOutput(args.out)
index = 0
if (args.operation == "enumerate"):
    comps = ABPUtils.GetComponents(graph,2)
    compIndex = 0
    
    for comp in comps:
        for c in comp:
            nameStr = ""
            groupStr = ""
            if ('name' in graph.node[c]):
                nameStr = "\t" + graph.node[c]['name']
            if ('color' in graph.node[c]):
                groupStr = "\t" + graph.node[c]['color']
                
            outFile.write(str(compIndex) + "\t" + str(c) + nameStr + groupStr + "\n")
        
        compIndex +=1
        




