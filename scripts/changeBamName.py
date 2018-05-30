#!/usr/bin/env python
import argparse
import sys

parser = argparse.ArgumentParser(description="")
parser.add_argument("infile", nargs="?", help="input bam file",  type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument("outfile",nargs="?", help="output bam file", type=argparse.FileType('w'), default=sys.stdout)
parser.add_argument('-d', action="store_true", default=False)
args = parser.parse_args()
 
from pysam import AlignmentFile
import re


"@RG     ID:c96857fdd5   PU:pileup_reads SM:NO_CHIP_ID   PL:PACBIO       DS:READTYPE=SUBREAD;CHANGELISTID=2.3.0.0.140018;BINDINGKIT=100356300;SEQUENCINGKIT=100356200;FRAMERATEHZ=100;BASECALLERVERSION=2.3;InsertionQV=iq;DeletionQV=dq;SubstitutionQV=sq;MergeQV=mq;SubstitutionTag=st;DeletionTag=dt;Ipd=ip"

bam1 = AlignmentFile(args.infile)
header = bam1.header.to_dict()
print(header)
#header["RG"]= [ { "SM":"NO_CHIP_ID",  "PL":"PACBIO", "ID":"BLASR", "PU":"pileup_reads"} ]
outbam = AlignmentFile(args.outfile, "wb", header=header)# template = bam1)


for idx,read in enumerate(bam1.fetch(until_eof=True)):
	#readID = read1.query_name.split("/")
	#read1.query_name = readID[0] + "/" + str(idx) 
	read.tlen = read.reference_length
	#read.query_name = "pileup_reads" + "/" + str(idx)
	#read.tags += [("RG","NO_CHIP_ID")]
	outbam.write(read)

outbam.close()



