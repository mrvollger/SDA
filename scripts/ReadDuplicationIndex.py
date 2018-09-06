#!/usr/bin/env python
import argparse
parser = argparse.ArgumentParser(description="")
parser.add_argument('bams', nargs='+', help="input bam files to assign duplication indexs")
parser.add_argument('-o', '--out', default="read.duplication.index", help="duplication index file. Each row has three values, the read name, the duplication index, and then wheter the read could be uniquly phased to just one of the PSV clusters")
args = parser.parse_args()
import pysam
import re
from collections import Counter
#print(args.bams)

listofsets = []

def readBam(bamf):
	bam = pysam.AlignmentFile(bamf)
	names = set()
	for read in bam.fetch(until_eof=True):
		names.add(read.query_name)
	return(names)

for bamf in args.bams:
	listofsets.append( (readBam(bamf), bamf) ) 



dofsets = {}
allreads = []
idxs = []
for reads, filename in listofsets:
	match = re.match("group.(\d+)/H2.WH.bam", filename)
	if(match is None):
		continue
	index = int(match.group(1))
	idxs.append(index)
	dofsets[index] = reads
	allreads += list(reads)
idxs = sorted(idxs)
counts = Counter(allreads)

#print(counts)
#print( max(counts.values()))
out = ""
for idx in idxs:
	for read in dofsets[idx]:
		isunique = "NotUnique"
		if(counts[read] <= 1):
			isunique = "Unique"
		out += "{}\t{}\t{}\n".format(read, idx, isunique)

open(args.out, "w+").write(out)

