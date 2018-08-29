#!/usr/bin/env python
import sys
from Bio import SeqIO
from Bio import Seq

import argparse

<<<<<<< HEAD
ap = argparse.ArgumentParser(description="Merge contigs from sda run")
=======
ap = argparse.ArgumentParser(description="Merge contigs from abp run")
>>>>>>> c5f142477186657859a42b2bb5fa5f5213be97cf
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



    
