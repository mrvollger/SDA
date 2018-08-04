#!/usr/bin/env python
from Bio import SeqIO
from Bio import Seq
from Bio import SeqRecord

import sys
inFile = open(sys.argv[1])
outFile = open(sys.argv[2],'w')

seq = SeqIO.read(inFile,"fasta")
istart = 0
while (istart < len(seq)):
    if (seq[istart] != 'N'):
        break
    istart+=1
iend = len(seq)
while (iend > istart):
    if (seq[iend-1] != 'N'):
        break
    iend-=1
s2 = seq[istart:iend]

seq2 = SeqRecord.SeqRecord(s2.seq,id=seq.id, name="",description="")
SeqIO.write(seq2, outFile, "fasta")
