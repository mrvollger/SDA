#!/usr/bin/env python


import argparse
import ABPUtils
import numpy as np
import sys
ap = argparse.ArgumentParser(description="Sort by haplotype")
ap.add_argument("mat", help="snv matrix file")
ap.add_argument("pos", help="pos to sort on.", type=int)

args = ap.parse_args()

mat = open(args.mat)
freqLine = mat.readline()
freq = np.array(freqLine.split())

gtList = []
for line in mat:
    gtl = np.array(list(line.strip()))
    gtList.append(gtl)

#import pdb
#pdb.set_trace()
    
gt = np.array(gtList)

gti = np.argsort(gt[:,args.pos])

v = [''.join(gt[i]) for i in gti]

out = '\n'.join(v) + "\n"
fl = list(freqLine)
fl[args.pos] = '*'
sys.stdout.write(''.join(fl))
sys.stdout.write(out)



