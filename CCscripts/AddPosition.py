#!/usr/bin/env python

import argparse
import networkx as nx
import ABPUtils
import sys
import numpy as np
import itertools

ap = argparse.ArgumentParser(description="Prune non-solution edges from a graph")
ap.add_argument("--graph", help="Input graph", required=True)
ap.add_argument("--out", help="Output graph", required=True)
args = ap.parse_args()
g = ABPUtils.ReadGraph(args.graph)

pos = nx.spring_layout(g)


import pdb
pdb.set_trace()
