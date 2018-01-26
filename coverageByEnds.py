#!/usr/bin/env python
import argparse
import os
import sys
import re
import numpy as np
import intervaltree
import pandas as pd

parser = argparse.ArgumentParser(description="")
parser.add_argument("-a", "--reads", nargs="+", help="bed file(s) with read start and end locations" )
parser.add_argument("-b", "--regions", nargs= "+", help="bed file with regions to count within")
parser.add_argument("-o", "--out", help="bed file with numebr of starts and ends in each region")
args = parser.parse_args()

readFiles = args.reads
regionFiles = args.regions
regions = {}
df = None

def defineRegions():
	for myfile in regionFiles:
		f = open(myfile).readlines()
		for line in f:
			line = line.split()
			Chr = line[0]
			start = int(line[1])
			end = int( line[2] )
			if(Chr not in regions):
				regions[Chr] = intervaltree.IntervalTree()
			# first vlaue is number of starts in region, second is numer of ends in region
			regions[Chr][start:end+1] = [0,0]


def increment(point, Chr, startOrEnd):
	for region in regions[Chr][point]:
		region.data[startOrEnd] += 1


def addCounts(myfile):
	print(myfile)
	f = open(myfile).readlines()
	for line in f:
		line = line.split()
		Chr = line[0]
		if(Chr not in regions):
			continue
		start = int(line[1]) 
		end = int(line[2])
		increment(start, Chr, 0)
		increment(end, Chr, 1)

def readReads():
	for myfile in readFiles:
		addCounts(myfile)


def makeBed():
	global df
	out = ""
	for key in sorted(regions):
		tree = regions[key]
		for region in sorted(tree):
			out += "{}\t{}\t{}\t{}\t{}\t{}\n".format(key, region.begin, region.end-1, 
					region.data[0], region.data[1], max(region.data[0], region.data[1]) )
	open(args.out, "w+").write(out)
	df = pd.read_csv(args.out, header=None, sep="\t")
	df.columns = ["chr", "start", "end", "startCount", "endCount", "maxCount"]
	print( df[["startCount", "endCount", "maxCount"]].describe() )

def main():
	defineRegions()
	readReads()
	makeBed()

main()


sds = df.groupby(["chr"])["maxCount"].describe()

