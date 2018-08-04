#!/usr/bin/env python

import argparse
import random
import ABPUtils

ap = argparse.ArgumentParser(description="Run correlation clustering on a graph with + and - edges.")
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
ap.add_argument("--graph", help="Input graph", required=True)
ap.add_argument("--method", help="Method to run (default: cc-pivot)", default="cc-pivot")
args = ap.parse_args()


def CCPivot(g):
    pivot = g.node[random.ranint(0,len(g.node))]
    
    
    
