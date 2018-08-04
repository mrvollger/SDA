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

ABPUtils.WriteMETISFile(g, args.out)
