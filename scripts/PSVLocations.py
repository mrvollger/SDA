#!/usr/bin/env python
import argparse
import sys

parser = argparse.ArgumentParser(description="")
parser.add_argument("--bam", "-s", help="input bam file + index mapping contigs to location wehre original PSVs were called")
parser.add_argument("--psvs", "-p", nargs='+', help="List of vcf files describing the different PSVs, must have form, group.{\d+}.vcf",  type=argparse.FileType('r'))
parser.add_argument("--check", help="Add an input fasta file (ASM.assemblies.fasta) to check if the PSV values are correct", default = None)
parser.add_argument("outfile",nargs="?", help="output table file", type=argparse.FileType('w'), default=sys.stdout)
parser.add_argument('-d', action="store_true", default=False)
args = parser.parse_args()

import pysam
import re
import pandas as pd
import numpy as np


complement = {"A":"T", "T":"A", "G":"C", "C":"G", "N":"N"}

def readPsvs():
	dfs = []
	for myfile in args.psvs:
		match = re.match(".+\.(\d+)\.vcf", myfile.name)
		assert match is not None
		group = int(match.group(1))
			
		df = pd.read_table(myfile, comment="#", header = None)
		df = df[[0, 1, 3, 4]]
		df.columns = ["contig", "pos", "ref", "alt"]
		df["group"] = group
		# change from one to zero based indexing 
		df["pos"] = df["pos"] - 1
		df["group"] = group
		dfs.append(df)

	rtn = pd.concat(dfs)
	return(rtn)

def readAlns():
	samfile = pysam.AlignmentFile(args.bam)
	alns = {}
	for alnSeg in samfile.fetch(until_eof=True):
		qpos, rpos = zip(*alnSeg.get_aligned_pairs())
		rpos = np.array(rpos)
		qpos = np.array(qpos)
		alns[alnSeg.query_name]	= (alnSeg, rpos, qpos)
	return(alns)


def readPos(psvrow, algSeg):
	psvpos = psvrow["pos"]

	qpos = []; qval = []; qisrc = []; qname = []; qgroup = []	
	# get reads overallping at this exact pos
	for name in alnSeg:
		aln, rposes, qposes = alnSeg[name]
		mask = np.argmax(rposes == psvpos)
		rpos = rposes[mask]; pos = qposes[mask] # position of psv in query sequence 
		
		if(pos is not None):
			#print(rpos, pos, name, len(aln.seq))
			qname.append(name)
			match = re.match(".*group\.(\d+).*", name)
			assert match is not None
			qgroup.append( int(match.group(1)))
			
			if(aln.is_reverse):
				qisrc.append(True)
				qpos.append( len(aln.seq) - pos - 1)
				qval.append( aln.seq[pos].upper() ) # This is corrected to the complmented base at the end 
			else:
				qisrc.append(False)
				qpos.append(pos)
				qval.append(aln.seq[pos].upper() )



	df = pd.DataFrame({"qpos":qpos, "alt":qval, "ccid":qname, "group":qgroup, "isrc": qisrc})
	df["ref"] = psvrow["ref"]
	df["collapse"] = psvrow["contig"]
	
	# return none if 
	if(psvrow["group"] not in qgroup):
		return(None)

	# fix reverse complment
	truealt = []
	for idx, row in df.iterrows():
		if(row.isrc == True):
			truealt.append( complement[row.alt] )
		else:
			truealt.append(row.alt)
	df["truealt"] = truealt

	# get contigs that have the correct PSV and the correct CC group 
	group = df[ (df.group == psvrow["group"]) & (df.alt == psvrow["alt"]) ]	
	
	# ignore if PSv is found in multiple contigs (that are not in the same group)
	# numContigsWithPSV = group.shape[0]
	# numContigsWithAlt = (df.alt == psvrow["alt"]).sum()
	# if(numContigsWithAlt > numContigsWithPSV ):
	#	print("Multiple groups have the PSV in the assembly, ignoring PSV at {} for group {}".format(psvrow["pos"], psvrow["group"]), file=sys.stderr)
	#	return(None)
	# the above check seemes unessisary 
	
	group = group[["collapse", "ccid", "group", "qpos", "truealt"]]	
	return(group)
	
def checkPsvs(df):
	print("Checking if PSVs appear in assemblies", file=sys.stderr)
	from Bio import SeqIO
	recs = SeqIO.to_dict(SeqIO.parse(args.check, "fasta"))
	groups = df.groupby(by ="ccid")
	for name, group in groups:
		# skip if we cannot find fasta entry 
		if(name not in recs):
			continue 
		rec = recs[name]
		for idx, row in group.iterrows():
			pos = row["qpos"]
			alt = row["truealt"]
			recalt = rec.seq[pos].upper()
			#print(alt, recalt, pos, name, file=sys.stderr)
			assert alt == recalt, "PSV called inccorectly at {}:{}, {} instead of {}".format(name,pos, alt, recalt)



psvs = readPsvs() # dataframe of psvs
alnSeg = readAlns() # set of contigs mapped back to reference collapse (alignedSegment)
print("Completed reading of PSVs/Alignments", file = sys.stderr)
frames = []
for idx, row in psvs.iterrows():
	frame = readPos(row, alnSeg)
	if(frame is not None):
		frames.append(frame)

if(len(frames) > 0):
	df = pd.concat(frames)
	if(args.check is not None):
		checkPsvs(df)
	df.to_csv(args.outfile, header=None, index=False, sep="\t")
else:
	# if no results write empty file
	args.outfile.write("")





