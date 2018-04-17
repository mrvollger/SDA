#!/usr/bin/env python
import argparse
import os
import pandas as pd
import numpy as np
import sys

parser = argparse.ArgumentParser(description="Merge previously sorted bedfiles inot one bed file")
parser.add_argument('--merge', nargs='*', default=None)
parser.add_argument('--infile', "-i", nargs='?', type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument('--outfile', "-o", nargs='?', type=argparse.FileType('w'),default=sys.stdout)
args = parser.parse_args()
beds = args.merge
instream = args.infile
outstream = args.outfile


def readBed(filename):
	bed = pd.read_csv(filename, sep = "\t", header=None)
	names = []
	for i in range(len(bed.columns)):
		names.append(i)
	bed.columns = ["contig", "start", "end"] + names[3:] 
	return(bed)

def sortBed(df):
	df.sort_values(by=["contig", "start"], inplace=True)
	return 

def writeBed(df):
	df.to_csv(outstream, header=False, index=False, sep="\t", float_format='%.7f')


if(beds is not None):
	print(beds)
	for bed in beds:
		df = readBed(bed)
		break

else:
	df = readBed(instream)
	sortBed(df)
	writeBed(df)

