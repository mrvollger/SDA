#!/usr/bin/env python
import re
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("me", help="" )
parser.add_argument("mark",help="",default=None)
args = parser.parse_args()

me = open(args.me)
mark = open(args.mark)

def readToDict(f):
	myd = {}
	for line in f:
		line = line.strip()
		line = re.sub('/', '', line)
		line=line.split()
		if(line[0] not in myd):
			myd[line[0]] = []
		myd[line[0]].append("{}-{}".format(line[1], line[2]))

	for key in myd:
		myd[key] = sorted(myd[key])

	return(myd)


me = readToDict(me)
mark = readToDict(mark)


allkeys = set(me.keys() + mark.keys())

fme = open("mitchOnly.txt", "w+")
fmark = open("markOnly.txt", "w+")

for key in sorted(allkeys):
	s = "contig:\t" + key + "\n"
	if(key in me ):
		s += "mitch:\t" + "\t".join(me[key]) + "\n"
	else:
		s += "\n"
	
	if(key in mark):
		s += "mark:\t" + "\t".join(mark[key]) + "\n"
	else:
		s += "\n"

	if(key in me and key not in mark):
		fme.write(s + "\n")
	elif(key in mark and key not in me):
		fmark.write(s + "\n")
	print(s)


