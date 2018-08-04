#!/usr/bin/env python

import argparse
import sys

ap = argparse.ArgumentParser(description="Compare the 1st colum in a grouped read list to the last, the ground truth")
ap.add_argument("table", help="Group table")
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
args = ap.parse_args()

table = open(args.table)

clusters = {}
for line in table:
    v = line.split()
    idx = int(v[0])
    if (idx not in clusters):
        clusters[idx] = {}
    if (v[-1] not in clusters[idx]):
        clusters[idx][v[-1]]=0
    # count this now
    clusters[idx][v[-1]]+=1

# Now check group concordance

for c in sorted(clusters.keys()):
    total = 0

    for rg in clusters[c].keys():
        total += clusters[c][rg]
    if (total == 0):
        continue
    for rg in clusters[c].keys():
        print "{}\t{}\t{}\t{:2.2f}".format(c,rg,clusters[c][rg], float(clusters[c][rg])/total)

        
