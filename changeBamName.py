#!/bin/env python

from pysam import AlignmentFile
import re

bam1 = AlignmentFile("reads.bam")
outbam = AlignmentFile("reads.sample.bam", "wb", template = bam1)

for idx,read1 in enumerate(bam1.fetch(until_eof=True)):
    #readID = read1.query_name.split("/")
    #read1.query_name = readID[0] + "/" + str(idx) 
    read1.query_name = read1.query_name + "." + str(idx)
    outbam.write(read1)

outbam.close()

