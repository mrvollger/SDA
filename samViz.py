#/usr/bin/env python
import argparse

parser = argparse.ArgumentParser(description="vizualize snps in a sam file")
parser.add_argument("sam", help="sam file on which we want to display snps, either sam or bam" )
parser.add_argument("snp",help="vcf with snps in it",default=None)
parser.add_argument("--snp2", help="second snp file that will be colored differently", default=None )
parser.add_argument("--out", help="name of the output png, defualt is snpsInSam.png", default="snpsInSam.png" )
parser.add_argument('-d', action="store_true", default=False)
args = parser.parse_args()
SAM = args.sam
SNP = args.snp
SNP2= args.snp2
OUT = args.out
DEBUG=args.d

import glob
import subprocess
import os
import sys
import re
import itertools 
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from Bio import SeqIO
import pysam 
if sys.version_info[0] < 3: 
    from StringIO import StringIO
else:
    from io import StringIO

def runCmd(cmd):
    #print("starting", cmd)
    proc = subprocess.Popen( cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    err = err.decode("utf-8")
    out = out.decode("utf-8")
    if(len(err) > 1):
        print(cmd)
        print(err)
    #print("ending cmd", err, out)
    return(out)


class SamSnps():
    def __init__(self):
        self.name = None
        self.snps = None
        self.snps2 = None
        self.start = None
        self.end = None
    def __str__(self):
        return(self.name + ":" + str(self.snps))
   
# returns list of snp positions
def readVCF(filename):
    vcf = open(filename)
    snps = {}
    #CHROM  POS ID  REF ALT
    for line in vcf:
        if(line[0] == "#"):
            continue
        token = line.split("\t")
        # convert from 1 based to 0 based system 
        pos = int(token[1]) - 1
        ref = token[0]
        alt = token[4]
        #entry = "{}_{}_{}".format(token[0], pos, token[4])
        if(ref not in snps):
            snps[ref]={}
        snps[ref][pos]=alt
    if(DEBUG):
        for key in snps:
            print(sorted(snps[key]), len(snps[key]))
    return(snps)

def samLine(row, snps):
    ref = row.reference_name
    refSeq=row.get_reference_sequence()
    rstart = row.reference_start
    rend = row.reference_end
    rtn = [] 
    for pos in sorted(snps[ref]):
        if(rstart > pos or rend <= pos):
            continue
        rpos = pos - rstart
        alt = snps[ref][pos]
        if(refSeq[rpos] == alt):
            #print(pos, alt, "ref")
            rtn.append(pos)
    return(rtn)

def readSAM(samfile, snps):
    samfile = pysam.AlignmentFile(samfile)
    samsnps = []
    for idx, row in enumerate(samfile):
        read = SamSnps()
        read.snps = samLine(row, snps)
        read.name = row.query_name
        read.start = row.reference_start
        read.end = row.reference_end
        samsnps.append(read)
    if(DEBUG):   
        for read in samsnps:
            print(read)
    return(samsnps)


def addsnps(snps, ax, idx):
    for snp in snps:
        width = 50
        p = patches.Rectangle((snp, idx),   # (x,y)
            width,          # width
            1,          # height
            facecolor="red", linewidth=.25
            )
        ax.add_patch(p)
   

def plotSamSnps(read, ax, idx):
    # resize for a longer read
    if(read.end + 1000 > ax.get_xlim()[1] ):
        ax.set_xlim( [0, read.end+1000] )

    # add rectanlge for read 
    ax.add_patch(patches.Rectangle(
            (read.start, idx),   # (x,y)
            read.end - read.start,          # width
            1,          # height
        ))
    
    # add patches for snps
    addsnps(read.snps, ax, idx) 

def plotSnps(samsnps):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.set_ylim([0, len(samsnps)])
    ax.set_ylabel('Sam Read Index')
    ax.set_xlabel('Reference Position')

    for idx, read in enumerate(samsnps):
        plotSamSnps(read, ax, idx)
    fig.savefig(OUT, dpi=300, bbox_inches='tight')

def main():
    snps = readVCF(SNP)
    samsnps  = readSAM(SAM, snps)
    plotSnps(samsnps)


main()

