#!/bin/env python 
import glob
import subprocess
import argparse
import os
import re
from Bio import SeqIO
#os.system("module load pandas/0.13.1")
import numpy as np
import pandas as pd
import sys
if sys.version_info[0] < 3: 
    from StringIO import StringIO
else:
    from io import StringIO


parser = argparse.ArgumentParser()
parser.add_argument("--reads", help="reads file fasta, not used anymore ")
parser.add_argument("--ref", help="reference that I want to mask repeats from")
parser.add_argument("--alldir", help="update the contig files for all of the bacs", action="store_true")
parser.add_argument("--mapQ", help="quality threshold", default="0")
parser.add_argument("--maxM", help="try to force self mapping", default="20")
parser.add_argument("--nproc", help="number of cores to use, defualt one", default="1")
parser.add_argument("--maskSize", help="the min size a repeat must be to mask it", default="2000")
parser.add_argument("--out", help="output reference with masked repeats", default="changeme")
args = parser.parse_args()

maxM = args.maxM
nproc = args.nproc
mapQ = args.mapQ
ref = args.ref
reads = args.reads # this is not used anymore 
alldir=args.alldir
minMaskSize = int(args.maskSize)
outfile = args.out 

pctAccCut = 95.0
eliminated = []


def runBlasr(refS, readS):
    ref = "tempRef.fasta"
    read = "tempRead.fasta"
    SeqIO.write(refS, ref, "fasta")
    SeqIO.write(readS, read, "fasta")

    blasr_options = " -maxMatch " + maxM + " -nproc " + nproc  + " -bestn 10 " + "-minMapQV " + mapQ + " -preserveReadTitle "
    blasr_command = "blasr " + read + " " + ref + blasr_options + "-header -m 4 "
    proc = subprocess.Popen( blasr_command , stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    TESTDATA=StringIO(out)
    df = pd.read_csv(TESTDATA, sep=" ")
    #print(df["pctsimilarity"])
    return df

def sortByLength(records):
    sort = []
    sortD = {}
    for rec in records:
        sortD[len(rec.seq)] = rec
    
    for key in sorted(sortD):
        #print(key)
        sort.append(sortD[key])
    print("length of sequences")
    for rec in sort:
        print(rec.id + ": " + str(len(rec.seq)))
   
    return sort


def maskSeqs(df, records,  eliminated):
    for elim in eliminated:
        df = df[df.tname != elim]
    df = df[df.pctsimilarity > pctAccCut]
    df["maplength"] = df.tend - df.tstart
    df = df[df.maplength > minMaskSize]
    if(not df.empty):
        print("Alignments that call for masking: ")
        print(df.to_string())

    for index, row in df.iterrows():
        rec = records[0]
        counterCheck = 0
        for x in records:
            if(x.id == row.tname):
                rec = x
                counterCheck += 1
        assert(counterCheck == 1)
        
        replaceSeq = rec.seq[0:row.tstart] + "N" * row.maplength + rec.seq[row.tend:]
        if(row.tstrand == 1):
            end = row.tseqlength - row.tstart
            start = row.tseqlength - row.tend
            replaceSeq = rec.seq[0:start] + "N" * row.maplength + rec.seq[end:] 
        rec.seq = replaceSeq
    

def mask(records,  contigIdx):
    global eliminated
   #  and eliminate the first record from being comapred against
    smallRec = records[contigIdx]
    eliminated.append(smallRec.id)
    df = runBlasr(records, records[contigIdx])
    maskSeqs(df, records, eliminated)
    #print(eliminated)

def runForRef(ref):
    global outfile 
    records = list(SeqIO.parse(ref, "fasta"))
    # order the recrods 
    records = sortByLength(records)
    for idx, rec in enumerate(records): 
        mask(records, idx)
    if(outfile=="changeme"):
        outfile = "/".join(ref.split("/")[:-1] ) + "/" + "p_ctg.masked.fa"
    print("Output file: " + outfile)
    SeqIO.write(records, outfile, "fasta")


if(alldir):
    bacglob = "bacs/C*/unitigs/"
    for mydir in glob.glob(bacglob):
        thisRef = mydir + "p_ctg.fa"
        print("")
        print("Current Reference: " + thisRef)
        runForRef(thisRef) 
        eliminated = []
    #runForBac(outdir)
elif(ref):
    runForRef(ref)

os.system("rm tempRead.fasta")
os.system("rm tempRef.fasta")


