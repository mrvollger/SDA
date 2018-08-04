#!/usr/bin/env python


import argparse
import ABPUtils
import numpy as np
import sys
import networkx as nx
import matplotlib.pyplot as plt
import pickle


ap = argparse.ArgumentParser(description="Convert a weighted graph to a METIS file.")
ap.add_argument("graph", help="Input graph file.")
ap.add_argument("--out", help="Output file", default="/dev/stdout")

args = ap.parse_args()

g = ABPUtils.ReadGraph(args.graph)
nx.set_edge_attributes(g, 'weight', { (e[0],e[1]):-1*e[2]['capacity'] for e in g.edges(data=True)} )

mst = nx.minimum_spanning_tree(g)
nx.set_edge_attributes(mst, 'weight', { (e[0],e[1]):-1*e[2]['weight'] for e in mst.edges(data=True)} )    

ABPUtils.WriteGraph(mst, args.out)
