#!/usr/bin/env python

import argparse
import networkx as nx
import ABPUtils
import sys

ap = argparse.ArgumentParser(description="Converg sg_edges_list to graph")
ap.add_argument("input", help="Input file (should be sg_edges_list, but not assuming anything")
ap.add_argument("output", help="Output file. (gml format)")

args = ap.parse_args()


inFile = open(args.input)
if inFile is None:
    print "ERROR. Could not open " + args.input
    sys.exit(1)
edges = []
vertices = {}
for line in inFile:
    v = line.split()
    vertices[v[0]] = True
    vertices[v[1]] = True
    if (v[-1] != "TR"):
        edges.append((v[0],v[1]))

g = nx.Graph()
for v in vertices.keys():
    g.add_node(v)

for e in edges:
    g.add_edge(e[0],e[1])

ABPUtils.WriteGraph(g, args.output)
    
