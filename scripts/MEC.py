#!/usr/bin/env python
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("fragMat", help="assembly.counsensous.fragments.snv.mat" )
parser.add_argument("cuts", help="mi.gml.cuts" )
parser.add_argument("--out", help="output sam file with cc group parings", default = "reads.cc.bam" )
parser.add_argument("--sam", default="reads.bam", help="file with reas to be partitioned." )
parser.add_argument("-v","--verbose", action="store_true", default=False, help="verbose" )
parser.add_argument("-c","--compare", action="store_true", default=False, help="verbose" )
args = parser.parse_args()



import os
import sys
import re
import numpy as np
import glob
from collections import Counter
import pysam
if sys.version_info[0] < 3: 
    from StringIO import StringIO
else:
    from io import StringIO

m = 0
n = 0
NA = 2
psvPositions  = []

def readMat(myfile):
	mat = open(myfile)
	reads = {}
	global n, m  
	for read in mat.readlines():
		read = read.strip().split("\t")
		name = read[1]
		row = np.zeros(len(read[0]))
		for idx, pos in enumerate(read[0]):
			if(pos == "1" ):
				row[idx] = 1
			elif(pos == "."):
				row[idx] = 0
			elif(pos == "n"):
				row[idx] = np.nan
			else:
				print(pos)
		reads[name] = row	
	n = len(reads.keys())
	m = len(reads[reads.keys()[0]])
	return(reads)


def readCuts(myfile):
	cuts = open(myfile).readlines()
	haps = np.zeros((len(cuts), m))
	totalPSVs = 0
	for idx, psvs in enumerate(cuts):
		psvs=psvs.strip().split()
		psvs = map(int, psvs)
		for psv in psvs:
			psvPositions.append(psv)
		np.put(haps[idx], psvs, [1])		
		totalPSVs += len(psvs)
	# get rid of positions where there are two PSVs at one stop
	# not sure why this happens, will Have to fix with mark
	# confirmed as a bug, a fix is in the works
	# hack to skip it for now
	toRm = list(np.where( (haps.sum(axis = 0) > 1) )[0])
	for pos in toRm:
		print("mulitple PSVs in one spot!")
		psvPositions.remove(pos)

	return(haps)

def MEC(read, cuts, psvPos):
	hasValues = ~np.isnan(read)
	positions = psvPos & hasValues
	numPsvsPerCut = cuts[:, positions].sum(axis=1)
	minPossiblePSVs = min(numPsvsPerCut)
	readCount = read[positions].sum()
	readLength = positions.sum()
	resultMatrix = np.zeros([len(cuts), 7])
	# count the MEC for all groups
	for group, row in enumerate(cuts):
		# this is the MEC
		xor = np.logical_xor(read[positions], row[positions]).sum()
		# this represents the number of times both the cut and read had a shared PSV
		aand = np.logical_and(read[positions], row[positions]).sum()
		# maximum number of possible matches between read and cut (sum of the cut)
		cutCount = row[positions].sum()
		# number of PSVs 
		expectedMEC = readLength * 0.15
		resultMatrix[group, :] = [xor, aand, cutCount, group, readCount, readLength, expectedMEC] 	
	
	# require at least one match
	resultMatrix = resultMatrix[resultMatrix[:,1] > 0,:]
	out = ""
	group = -1 
	if(len(resultMatrix) > 0):
		# determine the best cut to assign the read to
		# get all of the lowest MEC values
		bestPoses = resultMatrix[:,0] == resultMatrix[:,0].min()
		# if there are ties get the one with the maximum number of matching PSVs
		bestPos = np.argmax(resultMatrix[bestPoses,1])
		# select the best group 
		best = resultMatrix[bestPoses,:][bestPos,:]	
		xor = best[0]; aand = best[1]; cutCount = best[2]; group = best[3]; eMEC = best[6]
		out += "MEC: {}/{}\tMatches: {}\tFracOfCut: {:.2f}\tGroup: {}\n".format(
					xor, eMEC, aand, aand/cutCount, group)
		out += str(resultMatrix) 

	elif(minPossiblePSVs == 0):
		out += "no possible assignment\n"
	else:
		out += "read possible originates from the original collapse\n"
	
	if(args.verbose):
		print(out + "\n")
	return(int(group))

def runMEC(reads, cuts):
	# areas where psvs actaully exist
	psvPos = np.zeros(m)
	psvPos[psvPositions] = 1
	psvPos = psvPos.astype(bool)
	# make sure that each column is unique to a cluster 
	#(make sure that there is one psv per column and only one)
	assert max(cuts[:, psvPos].sum(axis=0)) == 1
	assert min(cuts[:, psvPos].sum(axis=0)) == 1

	partition = {}	
	for idx, key in enumerate(sorted(reads)):
		read = reads[key]
		group = MEC(read,cuts, psvPos)
		if(group not in partition):
			partition[group] = []
		partition[group].append(key)
	
	#print(cuts[:,psvPos])
	return(partition)

def reverseDict(partition):
	readCuts = {}
	for group in partition:
		reads = partition[group]
		for read in reads:
			readCuts[read] = group
	return(readCuts)

np.set_printoptions(linewidth=3000)
reads = readMat(args.fragMat)
cuts = readCuts(args.cuts)
# key is the group number value is a list of reads for that group
partition = runMEC(reads, cuts)
# key is a read name, value is its group
readCuts = reverseDict(partition)
#del partition[-1]

inbam = pysam.AlignmentFile(args.sam)
outbam=pysam.Samfile(args.out, "wb", template=inbam)

# write sam files
for read in inbam.fetch(until_eof=True):
	name = read.query_name
	try:
		group = readCuts[name]
		read.tags += [('cc', group)]
		outbam.write(read)
	except Exception, err:
		print( name)
		print("There was no read by that name partitioned by MEC, most probably a bug")
		#exit()

inbam.close()
outbam.close()
exit()

# write to sam files
inbam = pysam.AlignmentFile(args.sam)
outsams = {}
invPar = {}
for group in partition:
	print(group, partition[group])
	name = "group." + str(group) + "/H2.MEC.bam"
	gsam = pysam.AlignmentFile( name , "w", template=inbam)
	outsams[group] = gsam
	# create opposite mapping of partition
	for read in partition[group]:
		invPar[read] = group

# write sam files
for read in inbam.fetch(until_eof=True):
	name = read.query_name
	if name in invPar:
		group = invPar[name]
		outsams[group].write(read)

if(args.compare):
	import pysam
	import re
	samfiles = sorted(glob.glob("group.*/H2.WH.bam"))
	sam = {}
	tnames = []
	for x in samfiles:
		ID = re.match("group.(\d+)/H2.WH.bam",x).group(1)
		ID = int(ID)
		samfile = pysam.AlignmentFile(x) 
		names = []
		for read in samfile.fetch(until_eof=True):
			name = re.match("(.*)\.\d+", read.query_name).group(1)
			names.append(name)
			tnames.append(name)
		sam[ID] = names
	tnames = Counter(tnames)
	

	out = "group\toverlap\tnumInH2\tnumInMEC\taveNumForDiff\taveForIsec\n"
	for key in sam.keys():
		x = set(sam[key])
		intersect = x.intersection(partition[key])
		diff = x.difference(partition[key])
		dtotal = 0.0
		itotal = 0.0
		for name in diff:
			dtotal += tnames[name]
		for name in intersect:
			itotal += tnames[name]
		dave = dtotal/(len(diff)+0.0001)
		iave = itotal/(len(intersect)+0.0001)
		
		out += "{}\t{}\t{}\t{}\t{:.2f}\t{:.2f}\n".format(
				key, len(intersect), len(sam[key]), len(partition[key]), dave, iave)
		
	print(out)	





