#!/usr/bin/env python
#!/usr/bin/env python
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("track",nargs='?', help="name of the track hub, CHM1, CHM13, etc", default="CHM1" )
args=parser.parse_args()
track = args.track

# define some things to start
import re
import os
import string 
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
				self.chr, self.start, self.end,
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


# write html
html = "<h2>{1} collpase from {0}</h2>".format(track, collapse) 
summary = open("summary.txt").read()
html += "<pre>" + summary + "</pre>"
html += "<a href=https://eichlerlab.gs.washington.edu/help/mvollger/{}/psvGraphs/{}.pdf>PSVgraph</a>".format(track, collapse)
html = string.replace(html, "\n", "<br>")
html = string.replace(html, "\t", " ")
html = getGenomeBrowserLinks(html)
#html = string.replace(html, " ", "&nbsp")
open(collapse + ".html", "w+").write(html)

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

m5 = open("WH_dup.m5")
for idx, line in enumerate(m5):
	token = line.strip().split(" ")
	key = token[5]
	temp = re.split(":|-", key)
	chr = temp[0]
	start = int(temp[1])
	end = int(temp[2])
	contig = token[0]
	mstart = int(token[7])
	mend = int(token[8])
	strand = token[9]
	diff = mend - mstart 
	perID = float(token[11])/(float(token[11]) + float(token[12]))
	realStart = start + mstart
	realEnd = realStart + diff

	rgb = "0,0,0"
	if(perID >= 0.998):
		rgb = "0,128,0"
	elif(perID >= 0.95):
		rgb="128,128,128"
	else:
		rgb="255,255,0"

	if(idx == 0):
		#print("browser position " + "{}:{}-{}".format(chr, start,end) )
		#print('track name=ItemRGBDemo type=bedDetail description="Item RGB demonstration" visibility=2 itemRgb="On"')
		foo = 1

	B = BedRec(chr=chr, start=realStart, end=realEnd, name=collapse, strand=strand, rgb=rgb, ID=contig, details=html)
	rtn += str(B) + "\n"


f=open("asm.bed", "w+")
f.write(rtn)
f.close()




#
# now from that information I want to make a bigbed file that will be acceptable as a track hub
#
#import pandas as pd
#names  = ["chr", "start", "end", "contig", "score", "strand", "x", "y", "rgb", "bc","bs", "be", "collapse", "description"]


# write bed file
#df = pd.read_csv("asm.bed", header=None, names=names, sep = "\t")
#bb = df[["chr", "start", "end", "contig", "score", "strand", "x", "y", "rgb"]]
#bb = df
#bb = bb.sort_values(['chr', 'start'], ascending=[True, True])
#bb.to_csv(collapse+".bedDetail", header=False, index=False, sep="\t")

# create big bed
#cmd = """fetchChromSizes hg38 > chrom.sizes && 
#cmd = """bedToBigBed {0}.bedDetail chrom.sizes {0}.bb""".format(collapse) 
#os.system(cmd)

# write new bed file
os.system("cp asm.bed " + collapse + ".bedDetail")

# write track hub
trackDB="""track {1}
type bigBed 9
shortLabel {1}
longLabel {1}
visibility full
bigDataUrl {0}/{1}.bb
html {0}/{1}.html
itemRgb on
parent {0}

""".format(track, collapse)
open(collapse + ".db", "w+").write(trackDB)




