#!/usr/bin/env python

import argparse
import ABPUtils
import networkx as nx

ap = argparse.ArgumentParser(description="")
ap.add_argument("--graph", help="Graph file.", required=True)
ap.add_argument("--minComponentSize", help="Only print components larger than this", type=int, default=10)
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
args = ap.parse_args()

g = ABPUtils.ReadGraph(args.graph)
components = [[c for c in comp] for comp in nx.connected_components(g)]

for comp in components:
	edgeColors = {}
	nodeColors = {}
	nEdges = 0
	if len(comp) < args.minComponentSize:
		continue
	for c in comp:
		if 'color' in g[c]:
			col = g[c]['color']
			if col not in nodeColors:
				nodeColors[col] = 0
			nodeColors[col] +=1	 
		for n in g[c]:
			if n < c:
				continue
			if 'color' in g[c][n]:
				col = g[c][n]['color']
				if col not in edgeColors:
					edgeColors[col] = 0
				edgeColors[col] +=1 
			nEdges +=1
	maxColor = 0
	if (len(edgeColors.values()) > 0):
		maxColor = max(edgeColors.values())
	maxNodeColor = 0
	if (len(nodeColors.values()) > 0):
		maxNodeColor = max(nodeColors.values())
	print "{}\t{:2.2f}\t{}\t{:2.2f}".format(len(comp), float(maxNodeColor)/ len(comp), nEdges, float(maxColor)/nEdges)
