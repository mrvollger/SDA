#!/usr/bin/env python
# make  simple vcf file used by haptree form a

# input:
# Line like this:
#chrom   pos     na      ng      nc      nt      nins    ndel
#000119F 34936   3       2       15      289     9       8       G,G,A,CC,A,A,A,C        .
# Output
#chr1	1000	snp1	A	G	.	PASS	*	GT:GQ	0/1/1/1:100


import argparse
from Bio import SeqIO
from Bio import Seq
from Bio import SeqRecord
ap = argparse.ArgumentParser(description="Transform a freq file to a vcf file.\n" +
                             "\tThe input frequency is output of MpileupToFreq.py\n"+
                             "\tchrom   pos     na      ng      nc      nt      nins    ndel\n"+
                             "\tThe output is a very simple vcf, without comments\n"+
                             "\tchr1	1000	snp1	A	G	.	PASS	*	GT:GQ	0/1/1/1:100\n")
ap.add_argument("--freq", help="Input file.", default="/dev/stdin")
ap.add_argument("--ref", help="Input ref file.", default="/dev/stdin")
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
args = ap.parse_args()

refFile = open(args.ref)
outFile = open(args.out, 'w')

ref = SeqIO.to_dict(SeqIO.parse(refFile, "fasta"))
nucIndex =  {'A':2,'C':3,'G':4,'T':5}
otherNucs = {'A':['C','G','T'], 'G':['A','C','T'], 'C':['A','G','T'],'T':['A','C','G']}

freq = open(args.freq)
idx = 1
for line in freq:
    v = line.split()
    if v[0] not in ref:
        continue
    pos = int(v[1])
    refPos = pos-1
    refNuc = ref[v[0]].seq[refPos].upper()
    otherMaxCount = 0
    otherMaxNuc = ''
    
    if refNuc not in ['A', 'C', 'G', 'T']:
        continue
    for on in otherNucs[refNuc]:
        c = int(v[nucIndex[on]])
        if (c > otherMaxCount):
            otherMaxCount = c
            otherMaxNuc = on

    outFile.write("\t".join(str(i) for i in [v[0], pos, "snp{}".format(idx), refNuc, otherMaxNuc, ".", "PASS", "*", "GT:GQ", "0/1/1/1:100"]) + "\n")
    idx+=1

    
    
    


