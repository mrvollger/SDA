#!/usr/bin/env python

import argparse
import ABPUtils
import networkx as nx
import numpy as np


ap = argparse.ArgumentParser(description="Generate the consensus within each componenet based on overlaps")
ap.add_argument("graph", help="Input graph")
ap.add_argument("mat", help="Distance matrix")
ap.add_argument("--minSupport", help="Minimum support for call.", type=int,default=10)
ap.add_argument("--minRatio", help="Minimum ratio of support", type=float, default=0.75)
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
ap.add_argument("--converge", help="Min changes to considered converge.", type=int,default=20)
args = ap.parse_args()

matFile = open(args.mat)
mat = ABPUtils.ReadGenotypeMatrix(matFile)
g = nx.read_gml(args.graph)

#IPython.embed()

comps = ABPUtils.GetComponents(g)

comp = comps[0]
m=mat['mat']

gap = 'n'
ref = '.'
alt = '1'
m=mat['mat']

def Update(m, comps, final=False):
    nChanged = 0    
    for comp in comps:
        for ri in comp:
            row = m[ri,]
            cols = np.where(row != gap)[0]
            overlaps = np.array(sorted(g[ri].keys()))
            if (len(overlaps) == 0):
                continue
            for col in cols:
                aln = m[overlaps,col]
                nr = len(np.where(aln == ref)[0])
                na = len(np.where(aln == alt)[0])
                prev = m[ri,col]
                if (nr + na < args.minSupport):
                    m[ri,col] = gap
                elif (float(na)/(na +nr) > args.minRatio):
                    m[ri,col] = alt
                elif (float(nr)/(na+nr) > args.minRatio):
                    m[ri,col] = ref
                else:
                    if (final):
                        m[ri,col] = gap
                if (prev != m[ri,col]):
                    nChanged+=1
    return (m,nChanged)

nChanged = args.converge
iter = 0
while  nChanged >= args.converge:
    (m,nChanged) = Update(m,comps)
    print "iter " + str(iter) + " " + str(nChanged)
    iter+=1

(m,nChanged) = Update(m,comps, True)




outFile = ABPUtils.OpenOutput(args.out)
ABPUtils.WriteGenotypeMatrix(mat, outFile)


        
    
    
