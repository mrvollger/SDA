#!/usr/bin/env python
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("--mi", help="this shoudl be mi.mi, this is the mutual information between all pairs", default="mi.mi" )
parser.add_argument("--gml", help="mi.gml, these are the things that actually turned into positive edges", default="mi.gml" )
parser.add_argument("--max", help="maximum number of shared read to be considered a repulsion edge", default=1 , type=int)
parser.add_argument("--out", help="repulsion edges", default="mi.repulsion")
parser.add_argument("--shared", help="this is the min number of shared for positive edge", default=5 , type=int)
parser.add_argument("--lrt", help="this and shared are not used to calcaulate possible repulsion edges, but edges that would have been repulsion edges without this script", default=1.5 , type=float)
parser.add_argument('-d', action="store_true", default=False)
args = parser.parse_args()
DEBUG=args.d

import numpy as np
import networkx as nx
from sets import Set

g = nx.read_gml(args.gml)
comps = nx.connected_components(g)


pos = nx.get_node_attributes(g,'pos')
compByPos = {}
for compNum, comp in enumerate(comps):
	for n in comp:
		compByPos[ int( g.node[n]["pos"] ) ] = compNum


#print(compByPos)

# get the min and max pos for each node
minMax = {}
for n in g.nodes():
	pos = int( g.node[n]["pos"] )
	cut = g.neighbors(n)
	
	poses = [ int( g.node[mynode]["pos"] )  for mynode in cut   ]
	minMax[pos] =  ( min(poses), max(poses) )



def possibleRepulsion(i, j):
	rtn = True
	if(i == j):
		return(False)
	elif(i not in compByPos ):
		return(False)
	elif(j not in compByPos ):
		return(False)
	elif(compByPos[i] != compByPos[j] ):
		return(False)
	elif( not ( (minMax[i][0] < j) & (j  < minMax[i][1] ) )  ):
		return(False)
	elif( not ( (minMax[j][0] < i) & (i  < minMax[j][1]) )  ):
		return(False)
	return(rtn)

# read though mi.mi
mi = open(args.mi)
repulsion = open(args.out, "w+")
rep = []
newrep = 0
oldrep = 1
for line in mi:
	line = line.split()
	i, j, lrt, ilrt, jlrt, iMinor, jMinor, shared = int(line[0]), int(line[1]), float(line[2]), float(line[3]), float(line[4]), int(line[5]), int(line[6]), int(line[7])

	if( possibleRepulsion(i,j) ):
		if(shared <= args.max):
			newrep += 1
			repulsion.write("{}\t{}\n".format(i,j))
			rep.append((i,j))
		if( shared <= args.shared or lrt <= args.lrt ):
			oldrep += 1
			#print(lrt, ilrt, jlrt, shared)
print(newrep, oldrep)


if(DEBUG):
	lines = open("mi.gml.sites").readlines()
	cuts = {}
	for cut, sites in enumerate(lines):
		sites = sites.split()
		for num in sites:
			cuts[int(num)] = cut

	out = {}
	for i,j in rep:
		if(i in cuts and j in cuts):
			a = cuts[i]
			b = cuts[j]
			l = [a,b]
			key = "{}_{}".format(min(l), max(l))
			if(key not in out):
				out[key] = 0
			out[key] += 1



	print(out)




