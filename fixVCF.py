#!/usr/bin/env python
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--vcf", help="vcf file")
parser.add_argument("--out", help="vcf file", default="out.vcf")
parser.add_argument("--bam", help="bam file", default="reads.bam")
args = parser.parse_args()

vcf = args.vcf
outvcf = args.out
bam = args.bam

import re
from pysam import AlignmentFile

f = open(vcf)
vcf = f.read()
f.close()
#vcf = re.sub("0\|1:100", "0|1:100-1", vcf)
vcf = re.sub("GT", "GT:PS", vcf)
vcf = re.sub("NA", "0", vcf)



# read in the right sample name 
outbam = AlignmentFile( bam )

# create a fake sample name 
sampleName = re.findall("SM:[^\s]+", outbam.text )
sampleName=list(set(sampleName))

assert len(sampleName)==1, "multiple sample names?"
sampleName = re.sub("SM:", "", sampleName[0])

# repalce the sample name
vcf = re.sub("sample",sampleName,vcf)



# make output vcf
f = open(outvcf, "w+")
f.write(vcf)


