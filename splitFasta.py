#!/usr/bin/env python
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("fasta", help="input fasta" )
parser.add_argument("dir",help="output dir",default=None)
parser.add_argument('--noHeader', action="store_true", default=False)
args = parser.parse_args()

from Bio import SeqIO


recs = list(SeqIO.parse(args.fasta, "fasta") )

for rec in recs:
	if(args.noHeader):
		rec.description = ""
		rec.name = rec.id
	SeqIO.write(rec, args.dir + "/" + rec.id, "fasta")



