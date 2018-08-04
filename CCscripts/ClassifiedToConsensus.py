#!/usr/bin/env python

import argparse
import ABPUtils
import numpy as np

import IPython
ap = argparse.ArgumentParser(description="Given an SNV mat that is categorized by the ground truth, print the consensus of each genotype")
ap.add_argument("mat", help="Input matrix")
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
ap.add_argument("--min", help="Minimum count", default=10, type=int)
ap.add_argument("--minFrac", help="Minimum fraction for consensus", type=float, default=0.7)

args = ap.parse_args()

matFile = open(args.mat)
genotypes = ABPUtils.ReadGenotypeMatrix(matFile)

mat = genotypes['mat']
g   = genotypes['groups']
gk  = g.keys()


def GetCounts(cd):
    if ('.' in cd):
        ref = cd['.']
    else:
      ref = 0
    if ('1' in cd):
        alt = cd['1']
    else:
        alt = 0
    if ('n' in cd):
        unk = cd['n']
    else:
        unk = 0
    return (ref,alt,unk)


for gki in gk:
    mg = mat[g[gki],]
    cons = ['A']*mg.shape[1]
    for i in range(0,mg.shape[1]):
        countmat = np.unique(mg[:,i],return_counts=True)
        counts = {countmat[0][i]: countmat[1][i] for i in range(0,len(countmat[0]))}
        (ref,alt,unk) = GetCounts(counts)
        if (ref+alt < args.min):
            cons[i] = 'n'
            continue
        rFrac = float(ref)/(ref+alt)
        aFrac = float(alt)/(ref+alt)
        if (rFrac > aFrac and rFrac > args.minFrac):
            cons[i] ='.'
        elif (aFrac > rFrac and aFrac > args.minFrac):
            cons[i] = '1'
        else:
            cons[i] = 'u'
    print ''.join(cons)

cm = np.zeros((mat.shape[1], 3))

for i in range(0,mat.shape[1]):
    countmat = np.unique(mat[:,i],return_counts=True)
    counts = {countmat[0][j]: countmat[1][j] for j in range(0,len(countmat[0]))}
    v = GetCounts(counts)
    cm[i][0] = v[0]
    cm[i][1] = v[1]
    cm[i][2] = v[2]

#IPython.embed()
    
    
    
