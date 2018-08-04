#!/usr/bin/env python

import argparse
import ABPUtils
import numpy as np
import sys
ap = argparse.ArgumentParser(description="Partition snvs into k-haplotypes")
ap.add_argument("snvs", help="SNV file, generated from FragmentToHapTree.py --frag file.frag --vcf file.vcf --out file.snv --fragments")
#ap.add_argument("vcf", help="VCF file, giving positions of all considered variants.")
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
ap.add_argument("--start", help="Restrict window to positions staring here.", default=None, type=int)
ap.add_argument("--end", help="Restrict window to positions ending here.", default=None, type=int)
ap.add_argument("--named", help="Write named output.", default=False, action='store_true')
ap.add_argument("--pos", help="Write the table of used SNVs here.", default=None)
ap.add_argument("--mat", help="Output matrix file.", required=True)
args = ap.parse_args()


snvFile = open(args.snvs)
#vcfFile = open(args.vcf)

    
snvs = [ABPUtils.ParseSNVLine(line) for line in snvFile ]
# Each row in snvPos are the positions for SNVs for read i
snvPos = [np.array(sorted(v[1].keys())) for v in snvs]
# Each row in snvVal are the values for SNVs for read i
snvVal = [np.array([v[1][i] for i in sorted(v[1].keys())]) for v in snvs]

# Store all positions that are used
vcfPosDict = {}
idx={'A':0,'C':1,'G':2,'T':3}


    
for i in range(0,len(snvPos)):
    for j in range(0,len(snvPos[i])):
        p = snvPos[i][j]
        if p not in vcfPosDict:
            vcfPosDict[p] = {}
        nuc = snvVal[i][j]
        if (nuc not in vcfPosDict[p]):
            vcfPosDict[p][nuc] = 0

        vcfPosDict[p][nuc] +=1

def GetRank(countDict):
    counts = sorted(countDict.values(),reverse=True)

    major,majorCount = 'N',0
    minor,minorCount = 'N',0

    if (len(countDict) > 0):
        for k,v in countDict.iteritems():
            if (v == counts[0]):
                major=k
                majorCount=counts[0]
                break
        if (len(countDict) > 1):
            for k,v in countDict.iteritems():
                if (v == counts[1]):
                    minor=k
                    minorCount=counts[1]
                    break
    return (major,majorCount,minor,minorCount)
            

# Store all positions
vcfPos = np.array(sorted(vcfPosDict.keys()))
if (args.start is not None):
    vcfPos = vcfPos[vcfPos > start]
if (args.end is not None):
    vcfPos = vcfPos[vcfPos <= end]

if (args.pos is not None):
    pos = open(args.pos, 'w')
    for p in vcfPos:
        (major,majorCount,minor,minorCount) = GetRank(vcfPosDict[p])
        pos.write("{}\t{}\t{}\t{}\t{}\n".format(p,major,majorCount,minor,minorCount))
    pos.close()

nPos = len(vcfPos)

genotypes = []
for i in range(0,len(snvPos)):
    genotype = np.array(['n']*nPos)
    idx = np.searchsorted(vcfPos, snvPos[i][(snvPos[i] < vcfPos[-1]) & (snvPos[i] >= vcfPos[0])])
    ridx = np.searchsorted(snvPos[i], vcfPos[idx])
    
    genotype[idx] = snvVal[i][ridx]
    genotype[genotype == '0'] = '.'
    genotypes.append(genotype)

genotypeStrs = []
for gi in range(0,len(genotypes)):
    if (args.named):
        genotypeStrs.append(''.join(genotypes[gi]) + "\t" + snvs[gi][0] )
    else:
        genotypeStrs.append(''.join(genotypes[gi]))

refCounts = []
varCounts = []
fracs = []
for i in range(0,nPos):
    col = [g[i] for g in genotypes]
    nRef = col.count('.')
    nAlt = col.count('1')
    if (nRef + nAlt > 0):
        fracRef = nRef / (float(nRef + nAlt))
        fracs.append(str(int(np.floor(fracRef*10))))
    else:
        fracs.append('n')
outFile = open(args.mat, 'w')
outFile.write('\n'.join(genotypeStrs)+ "\n")
