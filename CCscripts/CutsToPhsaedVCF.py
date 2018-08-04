#!/usr/bin/env python


import argparse
import ABPUtils
import numpy as np
import sys
import networkx as nx
#import matplotlib.pyplot as plt
import pickle


ap = argparse.ArgumentParser(description="Convert a weighted graph to a METIS file.")
ap.add_argument("cuts", help="The cuts")
ap.add_argument("pos", help="The position and status of each minor allele (ref or alt) for each snv in the phase graph.  The positions are 0 based, and so the VCF positions are pos+1")
ap.add_argument("vcf", help="Original vcf file.")

ap.add_argument("--base", help="Base for phased vcfs.", default="group")
ap.add_argument("--minComponent", help="Minimum cluster size", default=2, type=int)
ap.add_argument("--summary", help="Write summary of clusters here",default=None)
ap.add_argument("--ref", help="Reference name for vcf", default="ref")
ap.add_argument("--reflen", help="Lenght of reference for vcf", default=1,type=int)
ap.add_argument("--subgraph", help="Write the subgraph here.", default=None)
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

vcf = {}
header = []
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
    vcf[int(v[1])] = [v[2],v[3],v[4],v[5],v[6],v[7],v[8]]

    

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
    phasedVcf.write("##contig=<ID={},length={}>\n".format(args.ref, args.reflen))
    phasedVcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n")

    for snv in sorted(list[comp]):
        snvPos = pos[snv][0]
        #
        # Is the haplotype the reference or alternative allele
        #
        hap = int(pos[snv][3])
        
        
        if (hap==0):
            genotype = "0|1"
        else:
            genotype = "1|0"

        if (snvPos + 1 not in vcf):
            sys.stderr.write("Source vcf missing entry " +
                             str(snvPos+1) +
                             ", however all entries should exist.\n")
            continue
        vcfPos = snvPos +1
        vcfEntry = vcf[vcfPos]
        infoStr="AN=2;AO={};RO={}AL={}".format(pos[snv][1],pos[snv][3],pos[snv][4])


        phasedVcf.write("{}\t{}\t.\t{}\t{}\t.\tPASS\t{}\tGT\t{}\n".format(contig, vcfPos, vcfEntry[1], vcfEntry[2],infoStr, genotype))
    phasedVcf.close()
    if (args.summary is not None):
        summary.write("{}\t{}\n".format(compIndex, len(comp)))



    if (args.subgraph is not None):
        cs = sorted(comp)
        sys.stderr.write(str(compIndex) + "\n")
        subgraphName = args.subgraph + "." + str(compIndex) + ".gml"
        subgraph = nx.subgraph(graph, cs)
        nx.set_node_attributes(subgraph, 'pos', dict(zip(cs, [pos[c][0] for c in cs])))
#n        import pdb
#        pdb.set_trace()
        ABPUtils.WriteGraph(subgraph, subgraphName)
        
    compIndex+=1
