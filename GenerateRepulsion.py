#!/usr/bin/env python
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("--mi", help="this shoudl be mi.mi, this is the mutual information between all pairs", default=None )
parser.add_argument("--gml", help="mi.gml, these are the things that actually turned into positive edges", default="" )
parser.add_argument("--max", help="maximum number of shared read to be considered a repulsion edge", default=1 , type=int)
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
	if(i not in compByPos ):
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
newrep = 0
oldrep = 1
for line in mi:
	line = line.split()
	i, j, lrt, ilrt, jlrt, iMinor, jMinor, shared = int(line[0]), int(line[1]), float(line[2]), float(line[3]), float(line[4]), int(line[5]), int(line[6]), int(line[7])

	if( possibleRepulsion(i,j) ):
		if(shared <= args.max):
			newrep += 1
		if( shared <= 5 or lrt <= 1.5):
			oldrep += 1
print(newrep, oldrep)

