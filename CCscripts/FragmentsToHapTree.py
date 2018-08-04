#!/usr/bin/env python

import argparse
import ABPUtils
import sys
ap = argparse.ArgumentParser(description="Transform a fragments file into the simple fragment tree format")
ap.add_argument("--frags", help="Fragment file.", default="/dev/stdin")
ap.add_argument("--vcf", help="VCF file.  Used to filter missing positions")
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
ap.add_argument("--fragment", help="Write the fragment name", action='store_true', default=False)
ap.add_argument("--pre", help="Minimum prefix-anchor",default=0,type=int)
ap.add_argument("--post", help="Minimum post-snv-anchor",default=0,type=int)
args = ap.parse_args()

outFile = open(args.out,'w')

def ParseSimpleVCFLine(vcfLine):
    v = vcfLine.split()
    return [v[0], int(v[1]), v[2], v[3], v[4], 0]

#
# Prepare vcf for filtering
# 
vcfFile = open(args.vcf)
vcf = [ParseSimpleVCFLine(line)for line in vcfFile]
vcfIndex = {}
for i in range(0,len(vcf)):
    vcf[i][-1] = i
    vcfIndex[vcf[i][1]] = i


vcfQuery = {v[1]:v for v in vcf}
#
# Transform fragments
#
print vcfIndex
fragFile = open(args.frags)
frags = [ABPUtils.ParseFragLine(line, True) for line in fragFile]
allSNVPos = {}
for f in frags:
    #{3: 0, 4: 1, 5: 1}
    pairs = []
    seen = {}
    for i in range(0,len(f[2])):
        snvPos = f[2][i]
        vcfPos = snvPos + 1
	if snvPos in seen:
	    continue
        seen[snvPos] = True
        if (vcfPos in vcfQuery):
            snvIndex = vcfQuery[vcfPos][-1]
            allele = '0'
            snvRead = f[5][i]
            snvAlt  = f[4][i]
            if (snvRead == snvAlt):
                allele = '1'
            if (f[6][i] >= args.pre and f[7][i] >= args.post):
		if snvPos + 1 in vcfIndex:
		    snvIndex = vcfIndex[snvPos+1]
                    pairs.append("{} : {}" .format(snvIndex, allele))
            allSNVPos[snvPos] = True
    if (len(pairs) > 0):
        fragName = ""
        if (args.fragment):
            fragName = f[0] + "\t"
        
        outFile.write(fragName + "{" + " , ".join(pairs) + "}\n")
sys.stderr.write(str(len(allSNVPos.keys())) + " snvs\n")
                    
