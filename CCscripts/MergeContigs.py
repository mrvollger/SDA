#!/usr/bin/env python
import sys
from Bio import SeqIO
from Bio import Seq

import argparse

ap = argparse.ArgumentParser(description="Merge contigs from abp run")
ap.add_argument("--assemblies", help="Assemblies", nargs="+")
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
args = ap.parse_args()

outFile = open(args.out, 'w')
for asm in args.assemblies:
    asmFile = open(asm)
    for seq in SeqIO.parse(asmFile,"fasta"):
        seq.id = asm
        seq.description=""
        seq.name=""
        SeqIO.write(seq, outFile, "fasta")



    
