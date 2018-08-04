#!/usr/bin/env python
import networkx as nx
import sys
import pickle
import argparse
import intervaltree
import bisect
import copy

import ABPUtils

ap = argparse.ArgumentParser(description="Make a graph from a .fragments file")
ap.add_argument("fragments", help=".fragmnes file: read n chrom [tuple#1... tuple#n], where tuple is in the format pos,ref,alt,read,pre,suf")
ap.add_argument("--overlaps", help="Write overlaps", default=None)
ap.add_argument("--minOverlap", help="Minimum fragment overlap", default=2,type=int)
args=ap.parse_args()
fragmentFile = open(args.fragments)

frags = np.asarray([ABPUtils.Fragment(line) for line in fragmentFile], dtype=object)

fragLen = np.asarray([f.Length() for f in frags])


frags2 = frags[fragLen > 0]

overlaps = ABPUtils.BuildOverlapGraph(frags2, 2)
