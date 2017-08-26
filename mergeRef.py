#!/bin/env python 
import argparse

parser = argparse.ArgumentParser(description="exmaple command: mergeRef.py --ref ref.fasta --out ref.merged.fasta")
parser.add_argument("--ref", help="a reference .fasta file")
parser.add_argument("--out", help=" output fasta file")
parser.add_argument("--pos", help=" output of contig start locaitons in new reference", default="contigStarts.txt")
args = parser.parse_args()

ref = open(args.ref)
outfasta = args.out


# set of a divider between contigs
divide = ""
divideLength = 5000
for i  in range(0,divideLength):
    divide += "N"

newRef = ">1F\n"
# skip first header line
ref.readline()
pos = len(newRef)
lastChars = []
for lineBad in ref:
    line = lineBad.strip()    
    if(line[0] == ">"):
        newRef += divide 
        pos += divideLength
        lastChars.append(pos)
        continue
    seqLen = len(line)
    #print(seqLen,pos)
    pos += seqLen
    newRef += line 

newRef += "\n"


out = open(outfasta, "w+")
out.write(newRef)
out.close()

contigStarts = ""
for pos in lastChars:
    contigStarts += str(pos) + "\n"
    #print(newRef[pos-1], newRef[pos])

f = open(args.pos, "w+")
f.write(contigStarts)
f.close()







