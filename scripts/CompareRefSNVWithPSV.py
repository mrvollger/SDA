#!/usr/bin/env python

import argparse
import numpy as np
from Bio import Seq
from Bio import SeqIO
from Bio import SeqRecord
ap = argparse.ArgumentParser(description="Compare PSVs derived from the reference to the psv graph")
ap.add_argument("--ref", help="Alignment of reference hits to the assembly")
ap.add_argument("--refFasta", help="Genome file.")
ap.add_argument("--vcf", help="Original vcf psv file.")
ap.add_argument("--psv", help="PSVs, one line per cluster.")
ap.add_argument("--nfti", help="read nucfreq.", default=None)
ap.add_argument("--nfto", help="Write nucfreq of matching positions.", default=None)
ap.add_argument("--writevcf", help="Output file.", default=None)

args = ap.parse_args()

ref = []
refFile = open(args.ref)
prevContig = None
posList = []
vcfFile = open(args.vcf)
vcf = []
vcfLines = []
refFastaFile = open(args.refFasta)
refSeq = SeqIO.read(refFastaFile, "fasta")

for line in vcfFile:
    if line[0] != "#":
        vcf.append(int(line.split()[1]))
        vcfLines.append(line.split())
        vcfLines[-1][7] = "AN=2;AO=NA;RO=NA;AL=NA"
        vcfLines[-1][8] = "GT:GQ"
        vcfLines[-1][9] = "0|1:100"

refCounts = {}
for line in refFile:
    vals = line.split()
    pos = int(vals[1])
    contig = vals[5]
    if pos + 1 not in refCounts:
        refCounts[pos+1] = 0
    refCounts[pos+1] +=1
    if contig != prevContig:
        if prevContig is not None:
            ref.append(posList)
        posList = []
    posList.append(pos+1)
    prevContig = contig
ref.append(posList)

multi = []
for rc in refCounts:
    if refCounts[rc] > 1:
        multi.append(rc)
multi.sort()
print len(multi)
for r in range(0,len(ref)):
    ref[r] = np.setdiff1d(ref[r], multi)

psvFile = open(args.psv)
psv = []

for line in psvFile:
    vals = line.split()
    psv.append([vcf[int(v)] for v in vals])

sims = []
for p in range(0,len(psv)):
    sim = [len(np.intersect1d(psv[p], ref[r])) for r in range(0,len(ref))]
    print str(p) + "\t" + " ".join(["{:3d}".format(d) for d in sim])
    sims.append(sim)

if args.nfti is not None or args.writevcf is not None:
    if args.nfti is not None:
        nfti = open(args.nfti)
        nfto = open(args.nfto,'w')
        nft = [ [int(v) for v in line.split()[1:]] for line in nfti ]
    
    for p in range(0,len(psv)):
        maxSim = max(sims[p])
        for i in range(0,len(sims[p])):
            if sims[p][i] == maxSim:
                break
        psvMatch = np.intersect1d(psv[p], ref[i])
        f = 0
        m = 0
        psvFilt = []
        nftFilt = []
        while (m < len(psvMatch) and f < len(vcfLines)):
            if psvMatch[m] == int(vcfLines[f][1]):
                psvFilt.append(vcfLines[f])
                if args.nfti is not None:
                    v = nftFilt[f] + [p]
                    nftFilt.append(v)

                m+=1
                f+=1
            elif psvMatch[m] < int(vcfLines[f][1]):
                m+=1
            else:
                f+=1
        if args.writevcf is not None:
            vcfOutName = args.writevcf + "." + str(p) + ".vcf"
            vcfOut = open(vcfOutName, 'w')
            vcfOut.write("""##fileformat=VCFv4.1
##INFO=<ID=AO,Number=A,Type=Integer,Description="Alternate allele observations">
##INFO=<ID=RO,Number=1,Type=Integer,Description="Reference allele observations">
##contig=<ID={},length={}>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample
""".format(refSeq.id, len(refSeq.seq)))
           
            vcfOut.write("\n".join(["\t".join(pf) for pf in psvFilt]) + "\n")
            vcfOut.close()
        if args.nfti is not None:
            nfto.write("\n".join(["\t".join(nf for nf in nftFilt)]) + "\n")                    
            

    if args.nfti is not None:
        nfto.close()

        
                
                
