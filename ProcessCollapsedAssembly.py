import os
import tempfile
import numpy as np
import pandas as pd
import json
import re
import glob
from pprint import pprint
from Bio import SeqIO

#
# setup the env for each exacution 
#
SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)
snake_dir = SNAKEMAKE_DIR + "/"
shell.executable("/bin/bash")
shell.prefix("source %s/env_python3.sh; set -eo pipefail; " % SNAKEMAKE_DIR)
RMenv = snake_dir + "env_RM.sh"

configFileName = "config/denovo.setup.config.json"
configfile: configFileName
reference = config["asm"]

RM_DB = "human"
if("species" in config): 
	RM_DB = config["species"].lower()
else:
	print('No "species" defined in the config for RepeatMasker, assuming {}.'.format(RM_DB) )


localrules: all, StartFofn, ConditionalOuts

report: "report/workflow.rst"


rule all:
	input:
		one="coverage/all.stats.txt",
		collapses="coverage/collapses.bed",
		json="config/sda.config.json",
		BedRgn="LocalAssemblies/BedAndRgn.done",
		sdarefs="LocalAssemblies/Refs.done",
		sdajson="LocalAssemblies/Configs.done",
		sdabams="LocalAssemblies/Bams.done",
		array="LocalAssemblies/RunAssembliesByArray.sh",
		extra = "LocalAssemblies/extra.done"
		#amb = os.path.dirname(os.path.realpath(reference)) + "/bwa.amb", # this is to index for bwa

DIRS = "reference fofns coverage LocalAssemblies alignments" 


#
# geting ready to run TRF and RM by splitting up the genome into multiple parts to run seperatly 
#
splitSize = 20 
FAI = config["asm"] + ".fai"
if(os.path.exists(FAI)):
	recs = open(FAI).readlines()
	if(splitSize > len(recs)):
		splitSize = len(recs)
else:
	recs = list(SeqIO.parse( config["asm"], "fasta"))
	if(splitSize > len(recs)):
		splitSize = len(recs)

#
# make an index for the masked assembly  
#
rule IndexASM:
	input:
		ref=reference,
	output:	
		fai=reference + ".fai",
	params:
		mem="8G",
		cores=1,
	shell:"""
samtools faidx {input.ref}
"""



rule splitRef:
	input:
		ref=reference,
		fai=ancient(FAI),
	output:
		split = temp( expand("reference/split/ref.{idx}.fasta", idx=range(0, splitSize) ) ),
		readme = "reference/README.txt",
	params:
		mem="16G",
		cores=1,
	run:
		shell("mkdir -p " + DIRS)
		shell("echo creating my own mask > " + output["readme"] )
		seqs = list(SeqIO.parse(input["ref"], "fasta"))
		toWrite = {}
		count = 0
		for idx, seq in enumerate(seqs):
			if(count not in toWrite):
				toWrite[count] = []
			toWrite[count].append(seq)
			count += 1
			if(count == splitSize):
				count = 0

		for key in toWrite:
			print(key, len(toWrite[key]))
			SeqIO.write(toWrite[key], output["split"][key], "fasta")
			# make a directory, becasue there are some weired thing with creating dirs
			shell("mkdir -p " + "reference/mask" + str(key))
#
#
#
rule RepeatMasker:
	input:
		split = "reference/split/ref.{idx}.fasta"
	output:
		RMout = temp("reference/mask{idx}/ref.{idx}.fasta.out"),
	params:
		mem="8G",
		cores=8,
	threads:8
	shell:"""
source {RMenv}
dir=reference/mask{wildcards.idx}
RepeatMasker \
	-species {RM_DB} \
	-e wublast \
	-dir $dir \
	-pa {threads} \
	{input.split}

"""

#
#
#
rule mergeRepeatMasker:
	input:
		split = expand("reference/mask{idx}/ref.{idx}.fasta.out", idx = range(0,splitSize))
	output:
		RMout = temp("reference/ref.RM.out"),
	params:
		mem="8G",
		cores=1,
	shell:
		"""
		cat {input.split} > {output.RMout}
		"""

#
#
#
rule RepeatMaskerBed:
	input:
		RMout = "reference/ref.RM.out",
	output:
		RM = temp("reference/ref.RM.out.bed"),
	params:
		mem="8G",
		cores=1,
	shell:
		"""
		{snake_dir}/scripts/RepeatMaskingToBed.py {input.RMout}
		"""
#
# merge the bed file from both TRF and RM programs 
# repeate makser along finds >98% of the repeats foudn by both RM and trf, 
# and I have found a bug in trf so I am droping it from the pipeline. 
#
rule mergeTRFandRM:
	input:
		#bed = "reference/trf.masking.bed",
		RM = "reference/ref.RM.out.bed",
	output:
		allR = "reference/all.repeats.bed",
	params:
		mem="8G",
		cores=1,
	shell:"""
cat {input} | sort -k1,1 -k2,2n | bedtools merge -i - > {output.allR}
"""

rule StartFofn:
	input: 
		reads = config["reads"], 
		fai=ancient(FAI),
	output:
		start = "fofns/startFofn"
	shell:
		"touch {output.start}"


nLines = len(open(config["reads"]).readlines())
maxIdx = int( np.ceil( nLines*1.0/config["read_files_per_job"] ) )
IDXs = list(range(maxIdx))
print("Using {} alignment jobs".format(maxIdx))

#
# split the baxh5 files into different fofns such that there are 10 bax files per fofn
#
rule SplitFOFN:
	input:
		fofn = config["reads"],
		mystart = "fofns/startFofn" #rules.StartFofn.output.start
	output:
		fofnSplit=expand("fofns/reads.{index}.fofn", index=IDXs)
	params:
		mem="1G",
		cores=1,
		ReadFilesPerJob=config["read_files_per_job"]
	run:
		fofn = open(config["reads"]).readlines()	
		nLines = len(fofn)
		batch = params["ReadFilesPerJob"]
		splitLists = [ fofn[i:(i + batch)] for i in range(0, nLines, batch) ]  
		for idx, sublist in enumerate(splitLists):
			out = open("fofns/reads.{}.fofn".format(idx), "w+")
			for line in sublist:
				line = os.path.abspath(line.strip()) + "\n"
				out.write(line)
			out.close()

#
#  For read depth, and other future steps, it is necessary to map reads back to the assembly.
#
MINALN=3000 # remove spurious alignments from common repeat elements 
if("minaln" in config):
	MINALN = config["minaln"]
print("Minimum Alignment Length: {}".format(MINALN))

MINSCORE=MINALN
ISPB=True; PBMM=False
if( "ont" in config ):
	if(config["ont"].lower() in  ["t", "true"]):
		ISPB=False

if( "pbmm2" in config ):
	if(config["pbmm2"].lower() in  ["t", "true"]):
		PBMM = True 
		DATATYPE = "SUBREAD"
	if("CCS" in config):
		if(config["CCS"].lower() in  ["t", "true"]):
			DATATYPE = "CCS"
		

if(ISPB and not PBMM):
	#
	# If a suffix arry does not yet exist for the assembly, build this.
	#
	rule MakeASMsa:
		input:
			asm=reference,
		output:
			asmsa=temp("reference/denovo.sa"),
		params:
			mem="16G",
			cores=1,
		shell:
			"sawriter {output.asmsa} {input.asm}"
	#
	# Map reads using blasr 
	#
	rule MapReads:
		input:
			asm=reference,
			asmsa="reference/denovo.sa",
			fofn="fofns/reads.{index}.fofn",
		output:
			align="alignments/align.{index}.bam"
		params:
			mem="8G",
			cores=8, 
		threads: 8
		shell:"""
blasr {input.fofn} {input.asm}  \
	--sa {input.asmsa} \
	--minAlnLength {MINALN} \
	--sdpTupleSize 13 \
	--maxMatch 25 \
	--bestn 1 \
	--advanceExactMatches 15  \
	--clipping subread \
	--nproc {threads} \
	--bam --out - | \
samtools view -b -F 4 - | \
samtools sort -@ {threads} -m 2G - -o {output.align}
"""

elif(PBMM):
	rule indexForPbmm2:
		input:
			reference,
		output:
			mmi=temp("reference/denovo.mmi"),
		params:
			mem="8G",
			cores=4,
		shell:"""
pbmm2 index --preset {DATATYPE} -j {params.cores} {input} {output}
"""
	
	rule MapReads:
		input:
			asm=reference,
			mmi="reference/denovo.mmi",
			fofn="fofns/reads.{index}.fofn",
		output:
			align="alignments/align.{index}.bam"
		params:
			mem="8G",
			cores=4, 
		threads: 8
		shell: """ 
pbmm2 align \
	--preset {DATATYPE} \
	-r 50000 \
	-j {threads} \
	{input.mmi} {input.fofn}  | \
	samtools view -bS -F 2308 - | \
	samtools sort -@ {threads} -m 2G -o {output.align} 
"""

else: # the reads are ont or not foramted like pb
	rule indexForMinimap:
		input:
			reference,
		output:
			mmi=temp("reference/denovo.mmi"),
		params:
			mem="16G",
			cores=1,
		shell:
			"""
			minimap2 -d {output} {input}
			"""

	rule MapReads:
		input:
			asm=reference,
			mmi="reference/denovo.mmi",
			fofn="fofns/reads.{index}.fofn",
		output:
			align="alignments/align.{index}.bam"
		params:
			mem="8G",
			cores=4, 
		threads: 8
		shell: """ 
minimap2 \
	-ax map-ont \
	-r 50000 \
	-s {MINSCORE} \
	--cs -L \
	-t {threads} \
	{input.mmi} $(cat {input.fofn}) | \
	samtools view -bS -F 2308 - | \
	samtools sort -@ {threads} -m 2G -o {output.align} 
"""


#
# index the alignments
# 
rule BamIndex:
	input:
		bam=rules.MapReads.output.align,
	output:
		bai="alignments/align.{index}.bam.bai"
	params:
		mem="4G",
		cores=2,
	threads: 2 
	shell:"""
samtools index -@ {threads} {input.bam}
"""

#
# make a bed version of each bam file, this will be used to calculate coverage
# No longer run, coverage is calcualted directly from the bam 
#
rule BamToBed:
	input:
		bam=rules.MapReads.output.align,
		bai=rules.BamIndex.output.bai,
	output:
		bed="alignments/align.{index}.bam.bed"
	params:
		mem="4G",
		cores=1,
	shell:"""
bedtools bamtobed -i {input.bam} | bedtools sort -i - > {output.bed}
"""

#
# this rule creats a bed file that is incremented by 100 for every contig
# these will be the feautes upon which we calculate depth wtih bedtools
#
rule FaiToBed:
	input:	
		asmfai=reference + ".fai",
	output:
		regions=temp("coverage/regions.100.bed"),
		regions1k=temp("coverage/regions.1000.bed"),
	params:
		mem="2G",
		cores=1,
	run:
		fai = open(input["asmfai"])
		out = ""
		outk = ""
		for line in fai:
			token = line.strip().split("\t")
			length = int(token[1])
			contig = token[0]
			cur = 0
			for nxt in range(100, length, 100):
				out += "{}\t{}\t{}\n".format(contig, cur, nxt-1)
				cur = nxt
			out += "{}\t{}\t{}\n".format(contig, cur, length-1)

			curk = 0
			for nxt in range(1000, length, 1000):
				outk += "{}\t{}\t{}\n".format(contig, curk, nxt-1)
				curk = nxt
			outk += "{}\t{}\t{}\n".format(contig, curk, length-1)

		outfile = open(output["regions"], "w+")
		outfile.write(out)
		open(output["regions1k"], "w+").write(outk)


#
# turn the bam files into a bed files that have the coverage of each 100bp segment of the genome
# 
rule BamToCoverage:
	input:
		bam=rules.MapReads.output.align,
		bai=rules.BamIndex.output.bai,
		regions=rules.FaiToBed.output.regions,
		#bed=rules.BamToBed.output.bed,
	output:
		#coverage="coverage/all.merged.bed",
		coverage=temp("coverage/coverage.{index}.bed"),
	params:
		mem="8G", # Tested on ~60X pacbio and it used >30G of ram, will use what is avalible
		cores=1,
	shell:"""
# get coverage and then sort by contig and then pos
bedtools coverage -bed -mean -sorted -a {input.regions} -b {input.bam} | \
		sort -k 1,1 -k2,2n > {output.coverage}
"""


#
# merge the coverage for each 100bp region in the genome 
#
rule MergeBed:
	input:
		coverage=expand("coverage/coverage.{index}.bed", index=IDXs),
	output:
		bedsort="coverage/all.merged.bed",
	params:
		mem="16G",
		cores=1,
	run:
		files = sorted(input["coverage"])
		cols = ["contig", "start", "end", "cov"]
		dtypes = {"contig":str, "start":int, "end":int, "cov":float}
		# set up shape of final dataframe
		shell("echo reading " + files[0])
		merged = pd.read_csv(files[0], sep="\t", header = None, names = cols, dtype=dtypes, engine="c" )
		
		for cur_file in files[1:]:
			shell("echo reading " + cur_file)
			df = pd.read_csv(cur_file, sep="\t", header = None, names = cols, dtype=dtypes, engine="c")
			merged["cov"] = merged["cov"] + df["cov"]
		
		merged.to_csv(output["bedsort"], sep="\t", header = False, index = False)


#
# calcualte the average coverages 
#
rule GetCoverageStats:
	input:
		one="coverage/all.merged.bed",
	output:
		one="coverage/all.stats.txt",
	params:
		mem='16G',
		cores=1,
	run:
		for infile, outfile in zip(input, output):
			bed = pd.read_csv(infile, sep = "\t", header=None,
					names=['contig', 'start', 'end',"coverage"])
			# I want to eliminte the really low or really high coverage things because they are probably
			# not assembled correctly and then assecess what the mean and standard deviation is
			top = bed.coverage.quantile(.95)
			bot = bed.coverage.quantile(.05)
			bed = bed[ ( bed.coverage < top ) & ( bed.coverage > bot ) ]
			
			stats = bed["coverage"].describe()
			out = "mean_coverage\tstd_coverage\n{}\t{}\n".format(stats["mean"], stats["std"])
			open(outfile,"w+").write(out)

#
# make configuration files that have the mincov, maxcov, and mintotal 
#
rule GenerateMinAndMax:
	input:
		stats="coverage/all.stats.txt"
	output:
		json="config/sda.config.json",
	params:
		mem='1G',
		cores=1,
	run:
		stats = pd.read_csv(input["stats"], header = 0, sep = "\t")
		# the plust one is to make it round up
		maxcov = int(stats.iloc[0]["mean_coverage"] )
		mintotal=int( stats.iloc[0]["mean_coverage"] + 3*stats.iloc[0]["std_coverage"] + 1 )
		# turns out sd is too varible to set a miniumum threshold.
		# mincov = int( stats.iloc[0]["mean_coverage"] - 3*stats.iloc[0]["std_coverage"] )
		mincov = int(stats.iloc[0]["mean_coverage"]/2.0)
		
		out2 = '{{\n\t"MINCOV" : {},\n\t"MAXCOV" : {},\n\t"MINTOTAL" : {},\n\t"project" : \"{}\"\n}}\n'.format(mincov, maxcov, mintotal, config["project"])
		open(output["json"], "w+").write(out2)

#
# count the number of ovlapping bases between the repeate masking and all,merged,.bed
#
rule CountOverlappingRepeatElements:
	input:
		combined="coverage/all.merged.bed",
		allR = "reference/all.repeats.bed",
	output:
		repeatCounted="coverage/all.repeatCounted.bed",
	params:
		mem='16G',
		cores=1,
	shell:
		"""
		bedtools intersect -a {input.combined} -b {input.allR} -wao | \
				 bedtools merge -i - -c 4,8 -o mean,sum > {output.repeatCounted}
				 # multiple rows will have the same region if there are two different repeate elements
				 # thus I need to get the mean of those two coverages (4,mean) and the the sum of the 
				 # overallping bases (8,sum)
		"""

def NoRepeatContent(row):
	val = "-"
	if(row["repeatCount"] <= 75):
		val = "+"
	return(val)

#
# confine bed file to only high coverage regions (hrc) and then merge any that are adj to eachother
#
rule BedForCollapses:
	input:
		combined="coverage/all.repeatCounted.bed",
		stats="coverage/all.stats.txt",
		fai=reference + ".fai",
	output:
		temp = temp("coverage/unmerged.tmp.collapses.bed"),
		collapses="coverage/unmerged.collapses.bed",
	params:
		mem='16G',
		cores=1,
	run:
		bed = pd.read_csv(input["combined"], sep = "\t", header=None, names=['contig', 'start', 'end', 'coverage', "repeatCount"])

		# change repete content to a percentage
		bed["repeatCount"] = 100 * bed["repeatCount"] / ( bed["end"] - bed["start"] )

		# require high enough coverage
		stats = pd.read_csv(input["stats"], header = 0, sep = "\t")
		#mincov= 2.0*stats["mean_coverage"] +3*stats["std_coverage"]
		mincov=int( stats.iloc[0]["mean_coverage"] + 3*stats.iloc[0]["std_coverage"] + 1 )
		print(mincov)
		bed = bed.ix[ bed["coverage"] >= mincov ]
		
		# marks the region as having or not having repeat content by strand
		bed["isNotRepeat"] = bed.apply(NoRepeatContent, axis=1)
		# writes to file before merging
		bed.to_csv(output["temp"], header=False, index=False, sep="\t" )
		# create a new file that merges adj regions that have the same strand (i.e. repeat content status) 
		shell("bedtools merge -i {output.temp} -d 2 -s -c 6,4,5 -o distinct,mean,mean > {output.collapses}")


#
# This function has high coverage regions made of repeate elements inhertet the coverage of adjacent
# unique high coverage regions
#
def removeIsolatedRepeatContent(df):
	rowNum = len(df)
	toKeep = [True]
	for idx in range(1, rowNum-1):
		notRepeat = list(df.iloc[idx-1:idx+2]["notRepeat"])
		keep = False
		if("+" in notRepeat ):
			keep=True 
		toKeep.append(keep)
	toKeep.append(True)

	df = df.ix[toKeep]
	return(df)
#
# This is a function for merge close high coverage regions 
#
def mergeHighCovRegions(df):
	rowNum = len(df)
	count = 0
	for idx in range(1, rowNum):
		row1 = df.iloc[idx-1]
		row2 = df.iloc[idx]
		tlength = row1["clength"] + row2["clength"]
		if( row1["contig"] == row2["contig"] ):
			start1 = row1["start"]; end1 = row1["end"]; start2 = row2["start"]; end2 = row2["end"]
			maxMergeDist = min( tlength/2 + 10, 10000 )
			if(start2 - end1 <= maxMergeDist ):
				df.iloc[idx-1] =  pd.Series({"contig":"remove","start":0, 
					"notRepeat":"-", "end":0, "coverage":0, "clength":0, "reapeatPer":0.0})
				notRepeat = "-"
				if(row1["notRepeat"] == "+" and row2["notRepeat"] == "+"): 
					notRepeat = "+"
				coverage = (row1["clength"]*row1["coverage"] + row2["clength"]*row2["coverage"])/(1.0*tlength)
				repeatPer= (row1["clength"]*row1["repeatPer"]+row2["clength"]*row2["repeatPer"])/(1.0*tlength)
				df.iloc[idx] =  pd.Series({"contig":row1["contig"], "start":start1, "end":end2, 
					"coverage":coverage, "clength":end2-start1+1, 
					"repeatPer":repeatPer, "notRepeat":notRepeat})
				count += 1
	newdf = df.ix[df["contig"] != "remove"] 
	return(newdf)

#
# this takes regions that are collapsed and puts them together if they are close enough to one another
#
rule MergeBedForCollapses:
	input:
		fai=reference + ".fai",
		collapses="coverage/unmerged.collapses.bed",
		allR = "reference/all.repeats.bed",
	output:
		unf = "coverage/unfiltered.collapses.bed",
	params:
		mem='8G',
		cores=1,
	run:	
		# read in the merged set
		HCR = pd.read_csv(input["collapses"], sep = "\t", header=None,
				names=['contig', 'start', 'end', "notRepeat", 'coverage', "repeatPer"])
		# calcualte collapse length, +1 is because they are inclusive ranges on both sides
		print(HCR.head())
		HCR["clength"] = HCR["end"] - HCR["start"] + 1
		# see the function description
		HCR = removeIsolatedRepeatContent(HCR)		
		
		# I think i should combine collapses that are within a certain distance of one another, maybe
		# this function does that, taking inot account the repeate content 
		print(len(HCR))		
		merged = mergeHighCovRegions(HCR)
		
		# read in the length of the contigs
		fai = pd.read_csv(input["fai"], sep = "\t", header=None, 
				names=['contig', 'length', 'OFFSET', 'LINEBASES', "LINEWIDTH"])
		fai = fai[["contig", "length"]]
		# this adds the contigs length to the collapse 
		merged = pd.merge(merged, fai, on='contig', how='inner')
		# creats a column that has the dist to either the beging or end of the contig, whichever is closer
		merged["distFromEnd"]=pd.concat([merged["start"], merged["length"]-merged["end"]], axis=1).min(axis=1)
		#merged = merged.ix[merged["distFromEnd"] <= 50000]
		
		# write unfiltered to file
		merged.to_csv(output["unf"], header=False, index=False, sep="\t" )



# the miniumum size of collapse in order to set it up for SDA
MINLEN = 15000
if("MinCollapseLength" in config):
	MINLEN = config["MinCollapseLength"]

rule FilterCollapses:
	input:
		unf = rules.MergeBedForCollapses.output.unf, 
		#unf = "unfiltered.collapses.bed",
	output:
		collapses="coverage/collapses.bed",
		#png ="coverage/SizeRepeatFilter.png",
	params:
		mem='8G',
		cores=1,
	run:
		collapses = pd.read_csv(input["unf"], sep = "\t", header=None, 
				names=['contig', 'start', 'end', "notR", 'coverage', "RC", "length", "contigl", "distToEnd"])
		minsize = MINLEN
		maxRC = 75
		
		# apply filter
		collapses = collapses.ix[(collapses["length"] >= minsize) & (collapses["RC"]<=maxRC)]
		outf = output["collapses"] 
		collapses.to_csv(outf, header=False, index=False, sep="\t" )

		# plot what filter will be 
		# removed R to make the install easier. So not generating this plot 
		#cmd = "{}scripts/PlotFilterBySizeAndRepeatContent.R --bed {} --png {} --size {} --repeatContent {}".format(
		#		snake_dir, input["unf"], output["png"], minsize, maxRC)
		#shell(cmd)








def GetRegionInfo(collapses):
	myfile = open(collapses)
	rtn = []
	for idx, line in enumerate(myfile):
		tokens = line.strip().split()
		chrm = tokens[0]
		start = tokens[1]
		end = tokens[2]
		#curdir = "LocalAssemblies/{}.{}.{}/".format(chrm, start, end)
		curdir = "LocalAssemblies/region_{}/".format(idx)
		rtn.append((curdir, chrm, start, end))
	myfile.close()
	return(rtn)

# creates a regions file that has all the regions that are collapsed
# and creates a directory for eahc one of those regions
#
rule LocalAssembliesRegions:
	input:
		collapses = rules.FilterCollapses.output.collapses,
		#collapses="collapses.bed",
	output:
		regions="LocalAssemblies/regions.txt",
	params:
		mem='8G',
		cores=1,
	run:
		info = GetRegionInfo(input["collapses"])
		regions = ""
		makedirs = ""
		for curdir, chrm, start, end in info:
			makedirs += "mkdir -p " + curdir + " && "
			match = re.match("LocalAssemblies/(region_\d+/)", curdir)
			assert match is not None, curdir
			region = match.group(1)
			regions += region + "\n"
		
		makedirs += 'echo "done making dirs"'
		shell(makedirs)
		open(output["regions"], "w+").write(regions)



#-------------------------------------------------------------------------------#
#
#  SETUP REQUIRED INPUT FILES FOR SDA
#
#-------------------------------------------------------------------------------#

def created(my_path):
	return( os.path.exists(my_path) and (os.path.getsize(my_path) > 0) )

#
# add a bed file to each region specifying where in asm they were 
#
rule LocalAssembliesBed:
	input:
		collapses="coverage/collapses.bed",
		regions="LocalAssemblies/regions.txt",
	output:
		done="LocalAssemblies/BedAndRgn.done",
	params:
		mem='8G',
		cores=1,
	run:
		info = GetRegionInfo(input["collapses"])
		for curdir, chrm, start, end in info:
			bed = "{}\t{}\t{}\n".format(chrm, start, end )
			rgn = "{}:{}-{}\n".format(chrm, start, end)
			open(curdir + "orig.bed", "w+").write(bed)
			open(curdir + "orig.rgn", "w+").write(rgn)
			assert created(curdir + "orig.bed")
			assert created(curdir + "orig.rgn")
		shell("touch {output.done}")



#
# using the .rgn files make a fasta file consisting of the collapse 
#
rule LocalAssembliesRef:
	input:
		done="LocalAssemblies/BedAndRgn.done",
		collapses="coverage/collapses.bed",
		asm=config["asm"]
	output:
		refs="LocalAssemblies/Refs.done",
	params:
		mem='4G',
		cores=1,
	run:
		info = GetRegionInfo(input["collapses"])
		for curdir, chrm, start, end in info:
			out = curdir + "ref.fasta"
			shell('samtools faidx {} "{}:{}-{}" > {}'.format(input["asm"], chrm, start, end, out))
			assert created(out)
		shell("touch {output.refs}")


# aggregate the output regions 
def getRegions():
	regions1, = glob_wildcards("LocalAssemblies/{region}/orig.rgn").region,
	regions2, = glob_wildcards("LocalAssemblies/{region}/orig.bed").region,
	assert len(regions1) == len(regions2)
	return(regions1)

def aggrBam(wildcards):
	regions = getRegions()
	tmp= expand("LocalAssemblies/{region}/reads.orig.bam", region = regions) 
	return(tmp)

#
# find the reads in all the bam files that map to that region
#
rule LocalAssembliesBam:
	input:
		done="LocalAssemblies/BedAndRgn.done",
		collapses="coverage/collapses.bed",
		alns=expand("alignments/align.{index}.bam", index = IDXs), 
	output:
		bams = "LocalAssemblies/Bams.done"
	params:
		mem='24G',
		cores=1,
	run:
		info = GetRegionInfo(input["collapses"])
		for curdir, chrm, start, end in info:
			out = curdir + "reads.orig.bam"
			prefix = curdir + "tmp." 
			for idx, bam in enumerate(input["alns"]):
				tmpbam = prefix + str(idx) + ".bam"
				shell('samtools view -b {} "{}:{}-{}" > {}; '.format(bam, chrm, start, end, tmpbam))			
				assert created(tmpbam)
			shell("samtools merge {} {}*.bam".format(out, prefix) )
			shell("rm {}*.bam".format(prefix) )
			assert created(out)

		shell("touch {output.bams}")

#
# copy over the min max stats so that SDA knows mincov, maxcov, mintotal
#
rule LocalAssembliesConfig:
	input:
		done="LocalAssemblies/BedAndRgn.done",
		json="config/sda.config.json",
		collapses="coverage/collapses.bed",
	output:
		#cov="LocalAssemblies/{region}/sda.config.json",
		covs="LocalAssemblies/Configs.done",
	params:
		mem='1G',
		cores=1,
	run:
		info = GetRegionInfo(input["collapses"])
		for curdir, chrm, start, end in info:
			out = curdir + "sda.config.json"
			shell("cp {} {}".format(input["json"], out))
			assert created(out)
		shell("touch {output.covs}")




arrayScript = """#!/usr/bin/env  bash
# this should be in the format -t 1-N_jobs
#$ -t 1-{}
# This is the number of concurrent jobs to run
#$ -tc 50
#
# The remainder are options passed to the script
#$ -S /bin/bash -V
#$ -P eichlerlab
#$ -l mfree=4G
#$ -l h_rt=06:00:00
#$ -pe serial 4
#$ -cwd
#$ -p -200
#$ -q eichler-short.q
#$ -o {}/LocalAssemblies/out/$JOB_NAME.$TASK_ID.o
#$ -e {}/LocalAssemblies/out/$JOB_NAME.$TASK_ID.e


cwd=`awk "NR == $SGE_TASK_ID"  regions.txt`
cd $cwd
{}/SDA
cd ..

"""
#
# This rule generates a script that submits and array job of SDAs
# it will not be helpful to those that do not have a sungird cluster
#
rule GenerateBatchRunScript:
	input:
		regions="LocalAssemblies/regions.txt",
	output:
		array = "LocalAssemblies/RunAssembliesByArray.sh",
	params:
		mem='1G',
		cores=1,
	run:
		path = os.getcwd()
		numJobs = len( open( input["regions"]).readlines() )	
		array =  arrayScript.format(numJobs, path, path, snake_dir) 
		open(output["array"], "w+").write(array)
		shell("mkdir -p LocalAssemblies/out; chmod 777 {output.array}")

#-------------------------------------------------------------------------------#










#-------------------------------------------------------------------------------#
#
#  SETUP OPTIONAL INPUT FILES FOR SDA
#
#-------------------------------------------------------------------------------#
if("reference" in config):
	#
	# create the sequences in the reference that best match what I am generating 
	#
	GRCh38 = config["reference"]
	fai = GRCh38 + ".fai"  
	sa = GRCh38 + ".sa"  
	#
	# combine all of the ref.fastas so I can align them with one command, (much faster)
	#
	rule combineRefFasta:
		input:
			collapses="coverage/collapses.bed",
			refs="LocalAssemblies/Refs.done",
		output:
			allref="LocalAssemblies/all.ref.fasta",
		params:
			mem='1G',
			cores=1,
		run:	
			info = GetRegionInfo(input["collapses"])
			shell("> {output.allref}")
			for curdir, chrm, start, end in info:
				fasta = curdir + "ref.fasta"
				shell("cat {} >> {}".format(fasta, output["allref"]))

	"""
	> {output.allref}
	for i in {input.dupref}; do
		echo $i
		cat $i >> {output.allref}
	done
	"""

	#
	# map the collapse the the human reference
	#
	rule duplicationsFasta:
		input:
			dupref=rules.combineRefFasta.output.allref,
		output:
			dupreffai="LocalAssemblies/all.ref.fasta.fai",
			dupsam="LocalAssemblies/all.ref.fasta.sam",
		params:
			mem='8G',
			cores=8,
		threads: 8
		shell:""" 
samtools faidx {input.dupref}
minimap2 \
	-ax asm20 \
	-N 30 -p .20 \
	--eqx \
	-r 50000 -s 10000 \
	-t {threads} \
	{GRCh38} {input.dupref} | \
	samtools view -h -F 4 - | \
	samtools sort -m 4G -T tmp -o {output.dupsam}
"""
		
	#
	# filter the sam file to only include high identity long contigs 
	#
	rule getHighIdentity:
		input:
			dupsam=rules.duplicationsFasta.output.dupsam,
			#dupsam="LocalAssemblies/all.ref.fasta.sam",
		output:
			duptsv="LocalAssemblies/all.ref.fasta.identity.tsv",
		params:
			mem='4G',
			cores=1,
		shell:
			"""
			{snake_dir}/scripts/samIdentity.py --header {input.dupsam} > {output.duptsv}
			"""
	#
	# generate two bed file for each assmebliy, one with the region of the collapse, and one with 100000 bp of slop 
	# on either side
	#
	rule ConvertTsvToBedAndRgn:
		input:
			duptsv=rules.getHighIdentity.output.duptsv,
		output:
			bedDone="LocalAssemblies/bed.done.txt",
		params:
			mem='4G',
			cores=1,
		run:
			#names=["contig", "start", "end", "read", "x", "y", "z", "z", "perID", "m", "mm", "i", "d"]
			#df = pd.read_csv( input["duptsv"], sep="\t", header=None, names=names)
			df = pd.read_csv( input["duptsv"], sep="\t")
			df["length"] = df["reference_end"] - df["reference_start"]
			df=df.ix[ (df["length"]>=5000) & (df["perID_by_events"] > 0.80) ]
			df.reset_index(drop=True, inplace=True)
			
			allbed = "LocalAssemblies/all.ref.fasta.bed"
			
			df[["reference_name", "reference_start", "reference_end"]].sort_values(by=['reference_name', 'reference_start']).to_csv(allbed, sep="\t", index=False, header=False)
			#shell("bedtools merge -d 5000 -i {} > {}".format(allbed + ".tmp", allbed ) )	
			
			shell("bedtools slop -i {} -g {} -b 100000 > {}".format(allbed, fai, allbed+".slop"))	
			slop = pd.read_csv( allbed + ".slop", sep="\t", header=None, names=["contig", "start", "end"])
			df["longstart"] = slop["start"].astype('int64')
			df["longend"] = slop["end"].astype('int64')

			grouped = df.groupby(["query_name"])
			counter = 0
			for name, group in grouped:
				counter += 1
				match =  re.search('(.+):(\d+)-(\d+).*', name)
				print(name, counter)
				if(match):
					region = "{}.{}.{}".format(match.group(1),match.group(2),match.group(3))
					# short bed file
					outfile = "LocalAssemblies/{}/ref.fasta.bed".format(region)
					group[["reference_name", "reference_start", "reference_end"]].drop_duplicates().to_csv(outfile, sep="\t", index=False, header=False)
					#shell( "sort -u {} -o {}".format(outfile, outfile) ) # remove duplicates 
					# long bed file
					outfile2 = "LocalAssemblies/{}/ref.fasta.long.bed".format(region)
					group[["reference_name", "longstart", "longend"]].drop_duplicates().to_csv(outfile2, sep="\t", index=False, header=False)
					#shell( "sort -u {} -o {}".format(outfile2, outfile2) ) # remove duplicates 
			shell("touch " + output["bedDone"])
		



	#
	# actaully fetch that region from the genome 
	#
	rule getReferenceSequences:
		input:
			bedDone="LocalAssemblies/bed.done.txt",
			regions="LocalAssemblies/regions.txt",
		output:
			refDone="LocalAssemblies/ref.done.txt",
		params:
			mem='4G',
			cores=1,
		run:
			regions = open(input["regions"])

			groups = []

			for region in regions.readlines(): 
				region = region.strip()[:-1]
				bedfile = "LocalAssemblies/{}/ref.fasta.long.bed".format(region)
				fastafile = "LocalAssemblies/{}/duplications.fasta".format(region)
				# the awk part gets rid of duplicates 
				cmd = "bedtools getfasta -fi {} -bed {} | awk '!a[$0]++' > {}".format(GRCh38, bedfile, fastafile)
				if(os.path.exists(bedfile)):
					groups.append(cmd)
			
			start = 0
			for end in range(0, len(groups), 50):
				cmd = " ; ".join(groups[start:end])
				shell(cmd)
				start = end
			cmd = " ; ".join(groups[start:len(groups)])
			shell(cmd)

			shell("touch " + output["refDone"])

	if("genes" in config):
		#
		# get intersection of genes
		#
		genes = config["genes"]
		rule intersectGenes:
			input:
				regions="LocalAssemblies/regions.txt",
				refdone = rules.getReferenceSequences.output.refDone,
			output:
				mydone = "LocalAssemblies/genes.done.txt"
			params:
				mem='4G',
				cores=1,
			run:
				regions = open(input["regions"])
				groups=[]
				for region in regions.readlines(): 
					region = region.strip()[:-1]
					bedfile = "LocalAssemblies/{}/ref.fasta.bed".format(region)
					outgenes = "LocalAssemblies/{}/ref.fasta.genes.bed".format(region)
					cmd = "bedtools intersect -wo -a {} -b {} | bedtools sort -i - > {} ".format(genes, bedfile, outgenes)
					if(os.path.exists(bedfile)):
						groups.append(cmd)
				
				start = 0
				for end in range(0, len(groups), 50):
					cmd = " ; ".join(groups[start:end])
					shell(cmd)
					start = end
				cmd = " ; ".join(groups[start:len(groups)])
				shell(cmd)
				
				shell("touch " + output["mydone"])

#-------------------------------------------------------------------------------#
















#############################################
### adding pilon/illumina to the pipeline ###
#############################################
'''
If you have two fastq files your command:
bwa mem -M -t 16 ref.fa read1.fq read2.fq > aln.sam
is absolutely fine. Your read2.fq is called mates.fq in the bwa examples. 
If you view the first lines of both files the read names are identical despite a 1 or 2 for the corresponding read of the pair.

If you only have one interleaved fastq file you would use the -p option:
bwa mem -M -t 16 -p ref.fa read.fq > aln.sam
In this case both reads of a pair are in the same fastq file successively. Have a look at the read names.
'''

if("illumina" in config):
	illumina = open( config["illumina"]).readlines()
	readidxs = set()
	for line in illumina:
		match = re.match(".*/([A-Za-z0-9]+)_[1|2].fastq.*", line)
		if(match):
			index = match.group(1)
			readidxs.add(index)
	readidxs = list(readidxs)
	print("Illumina read IDs:")
	print(readidxs)

	rule SplitIllumina:
		input:
			fofn = config["illumina"],
		output:
			read1=expand("illumina/reads/{index}_1.fastq.gz", index = readidxs),
			read2=expand("illumina/reads/{index}_2.fastq.gz", index = readidxs),
		params:
			mem='4G',
			cores=1,
		run:
			lines = open(input["fofn"]).readlines()
			for line in  lines:
				path = line.strip()
				outpath =  " " + os.getcwd() + "/illumina/reads/."
				shell("ln -s " + path + outpath)


	rule indexRefForBWA:
		input:
			asm=ancient(reference),
		output:
			amb = os.path.dirname(os.path.realpath(reference)) + "/bwa.amb",
			ann = os.path.dirname(os.path.realpath(reference)) + "/bwa.ann",
			bwt = os.path.dirname(os.path.realpath(reference)) + "/bwa.bwt",
			pac = os.path.dirname(os.path.realpath(reference)) + "/bwa.pac",
			sa = os.path.dirname(os.path.realpath(reference)) + "/bwa.sa",
		params:
			mem='24G',
			cores=1,
		run:
			asm_index = os.path.dirname(os.path.realpath(reference)) + "/bwa"	
			shell("bwa index -p " + asm_index + " -a bwtsw {input.asm}")


	#
	#  For read depth, and other future steps, it is necessary to map reads back to the assembly.
	#
	rule MapIllumina:
		input:
			read1=("illumina/reads/{index}_1.fastq.gz"),
			read2=("illumina/reads/{index}_2.fastq.gz"),
			amb = os.path.dirname(os.path.realpath(reference)) + "/bwa.amb",
			ann = os.path.dirname(os.path.realpath(reference)) + "/bwa.ann",
			bwt = os.path.dirname(os.path.realpath(reference)) + "/bwa.bwt",
			pac = os.path.dirname(os.path.realpath(reference)) + "/bwa.pac",
			sa = os.path.dirname(os.path.realpath(reference)) + "/bwa.sa",
		output:
			align="illumina/alignments/align.{index}.bam"
		params:
			mem='8G',
			cores=8,
		threads: 8
		run:
			asm_index = os.path.dirname(os.path.realpath(reference)) + "/bwa"	
			cmd = "bwa mem -M -t {threads} " + asm_index 
			cmd += " {input.read1} {input.read2} | samtools view -bS -F 4 - | "
			cmd += " samtools sort -@ {threads} -m 8G - -o {output.align}"
			shell(cmd)

	rule indexBWA:
		input:
			align=rules.MapIllumina.output.align
		output:	
			bai="illumina/alignments/align.{index}.bam.bai",
		params:
			mem='4G',
			cores=1,
		shell:
			"""
			samtools index {input.align}
			"""


	rule illuminaAlnDone:
		input:
			bai = expand(rules.indexBWA.output.bai, index = readidxs),
		output:
			illumina="illumina/aln.done.txt",
		params:
			mem='4G',
			cores=1,
		shell:
			"""
			touch {output.illumina}
			"""

	#
	# find the reads in all the bam files that map to that region
	#
	rule LocalAssembliesIlluminaBam:
		input:
			illumina=rules.illuminaAlnDone.output.illumina,
			refdone = rules.getReferenceSequences.output.refDone,
			bed = rules.LocalAssembliesRef.input.bed,
		output:
			bams="LocalAssemblies/{region}/illumina.orig.bam"
		params:
			mem='4G',
			cores=1,
		run:
			import pysam
			bed = open(input["bed"]).read().strip()
			token = bed.split()
			
			allreads = None
			for idx, bam in enumerate(sorted(glob.glob("illumina/alignments/align.*.bam"))):
				samfile = pysam.AlignmentFile(bam, "rb", check_sq=False)
				if(idx == 0 ):
					allreads = pysam.AlignmentFile(output["bams"], "wb", template=samfile)

				for read in samfile.fetch(token[0], int(token[1]), int(token[2])):
					allreads.write(read)
				samfile.close()		
			
			allreads.close()
			


	rule illuminaDone:
		input:
			local = dynamic(rules.LocalAssembliesIlluminaBam.output.bams),
		output:
			illumina="illumina/done.txt"
		shell:
			"""
			touch {output.illumina}
			"""






#
# CHECK FOR THE OPTIONAL OUTPUTS
#
def condInput(wildcards):
	illumina="illumina/done.txt"
	genes ="LocalAssemblies/genes.done.txt"
	dups="LocalAssemblies/ref.done.txt"
	rtn = []
	if("illumina" in config):
		rtn.append(illumina)
	if("genes" in config):
		rtn.append(genes)
	if("reference" in config):
		rtn.append(dups)
	return(rtn)

rule ConditionalOuts:
	input:
		cond = condInput,
	output:
		extra = temp("LocalAssemblies/extra.done"),
	params:
		mem='4G',
		cores=1,
	shell:"""
touch {output.extra}
"""
		




