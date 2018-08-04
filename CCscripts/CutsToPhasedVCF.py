#!/usr/bin/env python


import argparse
import ABPUtils
import numpy as np
import sys
import networkx as nx
#import matplotlib.pyplot as plt
import pickle


ap = argparse.ArgumentParser(description="Given cuts from MinDisagreeCut, make vcfs.")
ap.add_argument("cuts", help="The cuts")
ap.add_argument("pos", help="The position and status of each minor allele (ref or alt) for each snv in the phase graph.  The positions are 0 based, and so the VCF positions are pos+1")
ap.add_argument("vcf", help="Original vcf file.")

ap.add_argument("--base", help="Base for phased vcfs.", default="group")
ap.add_argument("--minComponent", help="Minimum cluster size", default=2, type=int)
ap.add_argument("--summary", help="Write summary of clusters here",default=None)
ap.add_argument("--ref", help="Reference FAI", default="ref")

args = ap.parse_args()

cuts = ABPUtils.ReadCuts(args.cuts)
posFile = open(args.pos)
vcfFile = open(args.vcf)
pos = []
for line in posFile:
    v = line.split()
    if (v[3] == 'N'):
        v[3] = 1
    pos.append([int(v[0]), int(v[1]), int(v[2]), int(v[3]), int(v[4])])

vcf = []
header = []
refFaiFile = open(args.ref)
refEntries = []
for line in refFaiFile:
    vals = line.split()
    refEntries.append([vals[0], vals[1]])

# lines makred by 05/10 were changed by mitchell, to fix the contig each PSV was being associated with
# start cahnge 05/10
contigNames = {}
# end change
lineNumber = 0
for line in vcfFile:
    if (line[0] == "#"):
        if (line.split()[0] != "#CHROM"):
            header.append(line)
            print header
        continue
    #
    # parse vcf lines
    #001936F 529     snp10   C       T       .       PASS    *       GT:GQ   0/1/1/1:100
    #
    v = line.split()
    # most important are the 2nd & 3rd entries
    contig = v[0]
    vcf.append([v[0], v[1], v[2],v[3],v[4],v[5],v[6],v[7],v[8]])
    # end change
    

# start change 05/10
print("unique contigs/chromosomes")
print(set(contigNames.values()))
# end change

compIndex = 0
if (args.summary is not None):
    summary = open(args.summary, 'w')
    summary.write("Group\tnComps\n")

for cut in cuts:

    outFileName = args.base + "." + str(compIndex) + ".vcf"
    phasedVcf = open(outFileName, 'w')
    
    phasedVcf.write("##fileformat=VCFv4.1\n")
    phasedVcf.write("##INFO=<ID=AO,Number=A,Type=Integer,Description=\"Alternate allele observations\">\n")
    phasedVcf.write("##INFO=<ID=RO,Number=1,Type=Integer,Description=\"Reference allele observations\">\n")
    # this isn't used in the next step so fake it for now.
    for r in refEntries:
        phasedVcf.write("##contig=<ID={},length={}>\n".format(r[0], r[1]))
    phasedVcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n")

    for vcfPos in sorted(list(cut)):

        #
        # Is the haplotype the reference or alternative allele
        #

        if (vcfPos >= len(vcf)):
            sys.stderr.write("Source vcf missing entry " +
                             str(snvPos+1) +
                             ", however all entries should exist.\n")
            continue
        vcfEntry = vcf[vcfPos]
        #infoStr="AN=2;AO={};RO={}AL={}".format(pos[vcfPos][1],pos[vcfPos][3],pos[vcfPos][4])
        # This isn't working right now, but these values aren't really useful anyway.
        infoStr="AN=2;AO=NA;RO=NA;AL=NA"
        genotype="0|1:100"
        # start cahnge 05/10
        #phasedVcf.write("{}\t{}\t.\t{}\t{}\t.\tPASS\t{}\tGT\t{}\n".format(contig, vcfPos, vcfEntry[1], vcfEntry[2],infoStr, genotype))
        phasedVcf.write("{}\t{}\t.\t{}\t{}\t.\tPASS\t{}\tGT\t{}\n".format(vcfEntry[0], vcfEntry[1], vcfEntry[3], vcfEntry[4],infoStr, genotype))
        # end chagne 

    phasedVcf.close()
    if (args.summary is not None):
        summary.write("{}\t{}\n".format(compIndex, len(cut)))



    compIndex+=1
