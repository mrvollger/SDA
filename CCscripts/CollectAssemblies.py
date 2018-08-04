#!/usr/bin/env python

from Bio import Seq
from Bio import SeqIO

import glob
import sys
import os

outFile = open(sys.argv[1], 'w')
name="consensus"
if len(sys.argv) > 2:
    name=sys.argv[2]

fileNames = glob.glob("*/group.*/assembly.{}.fasta".format(name))

for f in fileNames:
    if os.path.getsize(f) == 0:
        continue
    contigs = list(SeqIO.parse(open(f), "fasta"))
    if len(contigs) > 1:
	continue
    assembly = contigs[0]
    assembly.id=f
    assembly.title=""
    assembly.description=""
    SeqIO.write(assembly, outFile, "fasta")
    print f

outFile.close()
    
