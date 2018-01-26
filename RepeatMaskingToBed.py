#!/usr/bin/env python
import re
import argparse
import sys
if sys.version_info[0] < 3: 
	from StringIO import StringIO
else:
	from io import StringIO


parser = argparse.ArgumentParser(description="")
parser.add_argument("dat", help=".dat file from trf" )
parser.add_argument("--mat",help="output matrix file",default=None)
parser.add_argument("--bed",help="output bed file",default=None)
args = parser.parse_args()
dat = args.dat
bed = args.bed
mat = args.mat
if(bed is None):
	bed = dat + ".bed"
if(mat is None):
	mat = dat + ".mat"


def isUseLine(line):
	# check to see empty
	if(line==""):
		return("ignore")
	token = line.split()
	if( token[0].isdigit() ):
		return("region")
	elif( token[0] == "Sequence:"):
		return("seq")
	return("ignore")


def readDat(datfile):
	Seq=""
	outstr = ""
	dat = open(datfile)
	for line in dat:
		line = line.strip()
		status = isUseLine(line)
		if(status=="seq"):
			Seq = line.split()[1]
		elif(status == "region"):
			line = "\t".join( line.split() )
			outstr += "{}\t{}\n".format(Seq, line)
	return(outstr)


def toBed(mys):
	lines = mys.split("\n")
	rtn = ""
	for line in lines:
		line = line.strip()
		if(line == ""):
			continue 
		line = line.split("\t")
		rtn += "{}\t{}\t{}\n".format(line[0], line[1], line[2]) 
	return(rtn)

def RMtoBed(df):
	print("RM to bed file starting")
	b = df[["contig", "start", "end"]].copy()
	b.reset_index(drop=True, inplace=True)
	for idx in range(1, len(b)):
		e1 = b.at[idx-1, "end"]; c1 = b.at[idx-1, "contig"]
		s2 = b.at[idx, "start"]; c2 = b.at[idx, "contig"]
		if( c1 == c2 and e1 > s2 ):
			s1 = b.at[idx-1, "start"]; e2 = b.at[idx, "end"] 
			b.at[idx, "start"]  = s1
			b.at[idx, "end"] = max(e1, e2)
			b.at[idx-1, "contig"]  = "remove" 
			b.at[idx-1, "start"]  = 0
			b.at[idx-1, "end"]  = 0
		if( idx%100000 == 0):
			print("Done merging bed through line:\t{}".format(idx))
	b = b.ix[ b["contig"] != "remove" ]
	# write to a bed file	
	b.to_csv(bed, header=False, index=False, sep="\t")


def readRMout(dat):
	header = ["SWscore", "div", "del", "ins", "contig", "start", "end", "left",
			"indicator", "repeat", "class", "rstart", "rend", "rleft", "ID", "NotSure"]
	s = ""
	# remove the header and blank lines, and add a line so each rows has the smae number of columns 
	for idx, line in enumerate(open(dat).readlines()):
		line = "\t".join( line.strip().split() )
		if(len(line) == 0):
			continue
		if(line[0] not in ["S", "s"]):
			if(line[-1] == "*" ):
				s += line + "\n"
			else:
				s += line + " .\n"
		if( (idx+1)%1000000 == 0):
			print("Done preprocessing through:\t{}".format(idx))
	dat = StringIO(s)
	df = pd.read_csv(dat, header=None, names = header, sep="\t")
	
	print("sorting by contig and postion")
	df.sort_values(['contig', 'start'], inplace=True)
	RMtoBed(df)
	
	print("writing matrix file")
	df.to_csv(mat, index=False, sep="\t")





if(".dat" in dat):
	mys = readDat(dat)
	open(mat, "w+").write(mys)
	mybed = toBed(mys)
	open(bed, "w+").write(mybed)
elif(".out" in dat):
	import pandas as pd
	readRMout(dat)







