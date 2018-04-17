#!/usr/bin/env python
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("bam", help='regex of files to use, must be put in quotes, e.g. "test.*.sam" ' )
args = parser.parse_args()

import glob
import os
import sys
import re
from collections import Counter
import pysam 



def readSam(myfile):
	bam = pysam.AlignmentFile(myfile)
	nameLength = {}
	for idx,read in enumerate(bam.fetch(until_eof=True)):
		nameLength[read.query_name] = read.infer_query_length()
	return(nameLength)

names = []
reads = {}
for bam in glob.glob(args.bam):
	sys.stderr.write(bam + "\n")
	dic =  readSam(bam)
	reads.update(dic)
	names+=dic.keys()

counts = Counter(names)


bad = 0.0
badN = 0.0
good = 0.0
goodN = 0.0
for readid in counts:
	num = counts[readid]
	if(num > 1):
		badN += 1
		bad += reads[readid]
	else:
		goodN += 1
		good += reads[readid]

#print(goodN, good/goodN)
#print(badN, bad/badN)


print(badN/(goodN+badN))

