#!/usr/bin/env python


import argparse
import ABPUtils
import numpy as np
import sys
ap = argparse.ArgumentParser(description="Sort by haplotype")
ap.add_argument("mat", help="snv matrix file")
ap.add_argument("--out", help="Output file", default="/dev/stdout")
ap.add_argument("--alph", help="Alphabet to use, 3 characters: ref,alt,gap", default='.1n')

args = ap.parse_args()

alph     = list(args.alph)
mat      = open(args.mat)
outFile  = open(args.out, 'w')
#freqLine = mat.readline()
#freq     = np.array(freqLine.split())
#print freq
gtList = []
groups = {}
index = 0
groupList = []

for line in mat:
    v = line.split()
    gtl = np.array(list(v[0]))
    gtList.append(gtl)
    if (len(v) > 1):
        if (v[2] not in groups):
            groups[v[2]] = []
        groups[v[2]].append(index)
    index +=1
    groupList.append(v[2])
gt = np.array(gtList)

altList = []
refList = []

r=alph[0]
a=alph[1]
g=alph[2]
for i in range (0,gt.shape[0]):
    altList.append(np.where(gt[i] == a)[0])
    refList.append(np.where(gt[i] == r)[0])

ngt = len(gt)
# Compare distance to members in the group

allGroups = np.array(groups.keys())
def GetMean(a):
    if (len(a) == 0):
        return 0
    else:
        return np.mean(a)
for i in range(0,ngt):
    innerMat = []
    innerMis = []
    for j in groups[groupList[i]]:
        nMatch = len(np.intersect1d(altList[i],altList[j], assume_unique=True))
        nMis   = len(np.intersect1d(altList[i],refList[j], assume_unique=True))+\
         len(np.intersect1d(refList[i],altList[j], assume_unique=True))
        if (nMatch > 0 or nMis > 0):
            innerMat.append(nMatch)
            innerMis.append(nMis)
#        if (nMatch > 10 and nMis > 10):
#            print "="*len(gt[i])
#            print ''.join(gt[i])
#            print ''.join(gt[j])
#            print ""

    otherGroups = np.setdiff1d(allGroups, np.array([groupList[i]]))

    for g in otherGroups:
        outerMat = []
        outerMis = []
        for j in groups[g]:
#            print "i " + str(altList[i])
#            print "j " + str(altList[j])
            nMatch = len(np.intersect1d(altList[i],altList[j], assume_unique=True))
            nMis   = len(np.intersect1d(altList[i],refList[j], assume_unique=True))+\
              len(np.intersect1d(refList[i],altList[j], assume_unique=True))

            if (nMatch > 0 or nMis > 0):
                outerMat.append(nMatch)
                outerMis.append(nMis)
    
    outFile.write("{}\t{:2.2f}\t{:2.2f}\t{}\t{:2.2f}\t{:2.2f}\t{}\n".format(i, GetMean(innerMat), GetMean(innerMis), len(innerMat), GetMean(outerMat), GetMean(outerMis), len(outerMat)))
