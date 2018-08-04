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
ap.add_argument("--metis", help="Metis annotation.")
ap.add_argument("--out", help="Output file", default="/dev/stdout")



args = ap.parse_args()
metisFile = open(args.metis)

m = [ int(l.strip()) for l in metisFile ]

g = ABPUtils.ReadGraph(args.graph)

nx.set_node_attributes(g, 'metis', { i: m[i] for i in range(0,len(m)) })


ABPUtils.WriteGraph(g, args.out)
