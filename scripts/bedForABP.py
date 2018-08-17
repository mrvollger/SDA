#!/usr/bin/env python
import argparse
import pandas as pd
import re
import os
import string 

parser = argparse.ArgumentParser(description="")
parser.add_argument("--stats",  help="file with stats", default="abp.table.tsv" )
parser.add_argument("--summary", help="summary file")
parser.add_argument("--track", help="name of the track hub, CHM1, CHM13, etc", default="Mitchell_CHM1" )
parser.add_argument("--out", help="output file to write to" )
args=parser.parse_args()
track = args.track

# define some things to start
collapse = os.path.basename( os.getcwd())

class BedRec:
	def __init__(self, chr=None, start=None, end=None, 
			name="name", score="0", strand=".", 
			thickS=None, thinkE=None, rgb="0,0,0", ID=1, details="None"):
		self.chr = chr
		self.start = start
		self.end = end
		self.name = name
		self.score = score
		self.strand = strand
		if(thickS is None):
			self.thickS = self.start
			self.thickE = self.end
		else:
			self.thickS = thickS
			self.thickE = thickE
		self.rgb = rgb
		self.ID = ID
		self.details = details
		self.blockCount = 1
		self.blockStart = 0
		self.blockSizes = int(self.end) - int(self.start)
	
	def __str__(self):
		#rtn = (13*"{}\t"+"{}").format(
		rtn = (10*"{}\t"+"{}").format(
				self.chr, int(self.start), int(self.end),
				self.name, self.score, self.strand,
				self.thickS, self.thickE, self.rgb, 
				#self.blockCount, self.blockSizes, self.blockStart,
				self.ID, self.details)
		return(rtn)



def getGenomeBrowserLinks(mys):
	matches = re.findall("(chr.*?:[0-9]*-[0-9]*)", mys)
	linkFormat = "https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position={}%3A{}-{}&hgsid=626947147_3QGFLBYxBUzHUyZVc3x7Ivts3Dji"
	href= "<a href={}>{}</a>"
	for match in matches:
		pos = re.split(":|-", match)
		assert len(pos) == 3
		link = linkFormat.format(pos[0], pos[1], pos[2])
		link = href.format(link, match)
		mys = re.sub(match, link, mys) 
	return(mys)

#
# write html
#
html = "<h2>{1} collpase from {0}</h2>".format(track, collapse) 
summary = open(args.summary).read()
html += "<pre>" + summary + "</pre>"
html += "<a href=https://eichlerlab.gs.washington.edu/help/mvollger/{}/psvGraphs/{}.pdf>PSVgraph</a><br>".format(track, collapse)
html += "<a href=https://eichlerlab.gs.washington.edu/help/mvollger/{}/psvGraphs/{}.png>CoverageGraph</a><br>".format(track, collapse)
html = str.replace(html, "\n", "<br>")
html = str.replace(html, "\t", " ")
html = getGenomeBrowserLinks(html)
#open(collapse + ".html", "w+").write(html)



#
# read in reference regions
#
rtn = ""
refRegions = ""
chrs=[]
starts=[]
ends = []
rgbs = []
colapsePos = open("ref.fasta.bed") 
for line in colapsePos:
	token = line.strip().split("\t")
	chrs.append( token[0])
	starts.append(token[1])
	ends.append(token[2])
	rgbs.append("139,0,0")#dark red
	refRegions += "{}:{}-{},".format(token[0], token[1], token[2])

for idx, x in enumerate(starts):
	name = collapse.split(".")[0]
	rtn += str(BedRec(chr=chrs[idx], name=collapse, start=starts[idx],
		end=ends[idx], rgb=rgbs[idx],ID=name, details = html)) +"\n"



#
# read in assemblies 
#
df = pd.read_csv(args.stats, sep ="\t")
for index, row in df.iterrows():
	if(row["Status"] == "Failed"):
		continue 

	contig = row["query_name"]
	strand = "."
	perID = row["perID_by_matches"]
	chr = row["bestChr"]
	realStart = row["bestStart"]
	realEnd = row["bestEnd"]

	rgb = "0,0,0"
	grey=99.0
	print(collapse, perID)
	if(perID >= 100):
		rgb = "0,100,0"
	elif(perID >= 99.99):
		rgb = "0,150,0"
	elif(perID >= 99.8):
		rgb = "0,200,0"
	elif(perID >= grey):
		# covert perID to a fraction that I can multiply to get different levels of grey 
		frac = min( 100 - (perID-grey)/(99.8-grey), 1.0)
		value = int( 192*frac )
		rgb="{},{},{}".format(value, value, value)
	elif(perID >= 95.0):
		value = 230
		rgb="{},{},{}".format(value, value, value)
	else:
		rgb="255,255,0"

	B = BedRec(chr=chr, start=int(realStart), end=int(realEnd),
			name=collapse, strand=strand, rgb=rgb, ID=contig, details=html)
	rtn += str(B) + "\n"





f=open(args.out, "w+")
f.write(rtn)
f.close()

# this jsut copies the coverge plot so that I have a nice place to look at it in the html file
os.system("cp Coverage.png " + collapse + ".png")



