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
shell.executable("/bin/bash")
#shell.prefix("source %s/env_PSV.cfg; set -eo pipefail; " % SNAKEMAKE_DIR)
shell.prefix("source %s/env_PSV.cfg; " % SNAKEMAKE_DIR)

#
# A little complicated to find the temp dir
#
SSD_TMP_DIR = "/data/scratch/ssd"
if "TMPDIR" in os.environ:
    TMPDIR = os.environ['TMPDIR']
elif "TMPDIR" in config:
    TMPDIR = config['TMPDIR']
elif os.path.exists(SSD_TMP_DIR):
    TMPDIR = SSD_TMP_DIR
else:
    TMPDIR = tempfile.gettempdir()

configfile: "abp_on_denovo.json"

# if trf file exists mask the fasta file
#reference = "reference/ref.masked.fasta"
reference = config["asm"]
print("The reference is this file: {}".format(reference))


localrules: all, 
			#mask, no_mask,
			#MergeCoverage, 
			#FaiToBed, 
			#getCoverageStats, 
			#GenerateMinAndMax, 
			BedForCollapses,
			GetOneKregionCoverage,
			FiveKWindowStepOneK,
			MergeBedForCollapses, 
			ConvertTsvToBedAndRgn,
			#MergeBed, 
			#LocalAssembliesRegions,
			#LocalAssembliesBed,
			#LocalAssembliesRgn,
			#LocalAssembliesRef,
			#LocalAssembliesBam,
			#LocalAssembliesConfig,

rule all:
	input:
		#combined="coverage/all.merged.bed",
		#stats="coverage/all.stats.txt",
		#minmax="MinMax.sh",
		#collapses="collapses.bed",
		#readme = "reference/README.txt",
		#laregions="LocalAssemblies/regions.txt",
		#ref = "reference/ref.masked.fasta",
		#allsam="LocalAssemblies/all.ref.fasta.sam",
		#duptsv="LocalAssemblies/all.ref.fasta.identity.tsv",
		#cov=dynamic("LocalAssemblies/{region}/coverage.json"),
		refDone="LocalAssemblies/README.txt",

DIRS = "benchmarks reference fofns coverage LocalAssemblies alignments" 
#
# geting ready to run tfg by splitting up the genome into 10 parts to run seperatly 
#
rule splitRef:
	input:
		ref=config["asm"]
	output:
		split = expand("reference/split/ref.{idx}.fasta", idx=range(0,10) ),
		readme = "reference/README.txt",
	params:
		cluster=" -pe serial 1 -l mfree=16G -l h_rt=12:00:00",
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
			if(count == 10):
				count = 0

		for key in toWrite:
			print(key, len(toWrite[key]))
			SeqIO.write(toWrite[key], output["split"][key], "fasta")
			# make a directory, becasue there are some weired thing with creating dirs
			shell("mkdir -p " + "reference/mask" + str(key))
#
# Run trf on the splits of the genome
#
rule CreateMask:
	input:
		split = "reference/split/ref.{idx}.fasta"
	output:
		split = "reference/mask{idx}/ref.masked.{idx}.fasta"
	params:
		cluster=" -pe serial 1 -l mfree=16G -l h_rt=12:00:00",
	shell:
		"""
		dir=reference/mask{wildcards.idx}
		cd $dir
		trf ../../{input.split} 2 7 7 80 10 50 1999 -m -d
		mv *.mask ../../{output.split}
		cd ../../
		"""
#
#
#
rule RepeateMasker:
	input:
		split = "reference/split/ref.{idx}.fasta"
	output:
		RMout = "reference/mask{idx}/ref.{idx}.fasta.out"
	params:
		cluster=" -pe serial 4 -l mfree=4G -l h_rt=12:00:00",
	shell:
		"""
		dir=reference/mask{wildcards.idx}
		module load perl/5.14.2
		module load RepeatMasker/3.3.0
		RepeatMasker -e wublast \
				-species human \
				-dir $dir \
				-pa 4 \
				{input.split}
		"""

#
#
#
rule mergeRepeateMasker:
	input:
		split = expand("reference/mask{idx}/ref.{idx}.fasta.out", idx = range(0,10))
	output:
		RMout = "reference/ref.RM.out",
	params:
		cluster=" -pe serial 1 -l mfree=8G -l h_rt=12:00:00",
	shell:
		"""
		cat {input.split} > {output.RMout}
		"""

#
#
#
rule RepeateMaskerBed:
	input:
		RMout = "reference/ref.RM.out",
	output:
		RM = "reference/ref.RM.out.bed",
	params:
		cluster=" -pe serial 1 -l mfree=8G -l h_rt=12:00:00",
	shell:
		"""
		/net/eichler/vol2/home/mvollger/projects/abp/RepeatMaskingToBed.py {input.RMout}
		"""


#
# merge the maksed fasta files generated by trf
#
rule mergeMask:
	input:
		split = expand("reference/mask{idx}/ref.masked.{idx}.fasta", idx = range(0,10))
	output:
		ref = "reference/ref.masked.fasta",
	params:
		cluster=" -pe serial 1 -l mfree=8G -l h_rt=12:00:00",
	run:
		seqs = []
		for myfile in input["split"]:	
			seqs = seqs + list(SeqIO.parse(myfile, "fasta"))
		#seqs = sorted(seqs)
		SeqIO.write(seqs, output["ref"], "fasta")

#
# Make bed files
#
rule CreateBed:
	input:
		masked = "reference/mask{idx}/ref.masked.{idx}.fasta"
	output:
		bed = "reference/mask{idx}/ref.masked.{idx}.bed"
	params:
		cluster=" -pe serial 1 -l mfree=8G -l h_rt=12:00:00",
	shell:
		"""
		input=reference/mask{wildcards.idx}/*.dat
		/net/eichler/vol2/home/mvollger/projects/abp/RepeatMaskingToBed.py $input --bed {output.bed}
		"""
#
# merge the bed files generated by trf
#
rule mergeBed:
	input:
		bed = expand("reference/mask{idx}/ref.masked.{idx}.bed", idx = range(0,10))
	output:
		bed = "reference/trf.masking.bed",
	params:
		cluster=" -pe serial 1 -l mfree=8G -l h_rt=12:00:00",
	shell:
		"""
		cat {input.bed} > temp.bed
		sort -k1,1 -k2,2n temp.bed > temp.sorted.bed 
		bedtools merge -i temp.sorted.bed > {output.bed}
		rm temp.bed temp.sorted.bed
		"""

#
# merge the bed file from both repeate programs 
#
rule mergeTRFandRM:
	input:
		bed = "reference/trf.masking.bed",
		RM = "reference/ref.RM.out.bed",
	output:
		allR = "reference/all.repeats.bed",
	params:
		cluster=" -pe serial 1 -l mfree=8G -l h_rt=12:00:00",
	shell:
		"""
		cat {input.bed} {input.RM} > temp.bed
		sort -k1,1 -k2,2n temp.bed > temp.sorted.bed 
		bedtools merge -i temp.sorted.bed > {output.allR}
		rm temp.bed temp.sorted.bed
		"""


#
# If a suffix arry does not yet exist for the assembly, build this.
#
rule MakeASMsa:
    input:
        asm=reference
    output:
        asmsa=reference + ".sa"
    params:
        cluster=" -pe serial 1 -l mfree=48G -l h_rt=6:00:00",
        blasr=config["blasr"]
    shell:
        """
        {params.blasr}/alignment/bin/sawriter {input.asm}
        """
#
# make an index for the masked assembly  
#
rule IndexASM:
	input:
		ref=reference,
	output:	
		fai=reference + ".fai",
	params:
		cluster=" -pe serial 1 -l mfree=8G -l h_rt=1:00:00",
	shell:
		"""
		samtools faidx {input.ref}
		"""

#
# split the baxh5 files into different fofns such that there are 10 bax files per fofn
#
rule SplitFOFN:
    input:
        readme = "reference/README.txt",
        fofn=config["reads"]
    output:
        fofnSplit=dynamic("fofns/reads.{index}.fofn")
    params:
        cluster=" -pe serial 1 -l mfree=1G -l h_rt=1:00:00",
        baxPerJob=config["bax_per_job"]
    shell:
        """
        split --numeric --lines {params.baxPerJob} {input.fofn} reads.
        for f in `ls reads*`; do
        mv $f fofns/$f.fofn
        done
        """

#
#  For read depth, and other future steps, it is necessary to map reads back to the assembly.
#
rule MapReads:
    input:
        asm=reference,
        asmsa=reference + ".sa",
        fofn="fofns/reads.{index}.fofn"
    output:
        align="alignments/align.{index}.bam"
    params:
        cluster=" -pe serial 8 -l mfree=8G -l h_rt=24:00:00",
        blasr=config["blasr"]
    shell:
        """
		# bestn 2 is because a asm might be fragmented in half and we want
		# a read to be able to map to both halves.
		# not sure if mapQV should be 0 or 30, Mark said 0 was probably a msitake and to use 30
        {params.blasr}/alignment/bin/blasr {input.fofn} {input.asm}  -sa {input.asmsa} \
			-sdpTupleSize 13  -sdpMaxAnchorsPerPosition 10 -maxMatch 25 \
			-minMapQV 30 -bestn 2 -advanceExactMatches 15  \
            -clipping subread -nproc 8 -sam -out /dev/stdout | \
            samtools view -bS -F 4 - | \
            samtools sort -T {TMPDIR}/blasr -@ 8 -m 8G - -o {output.align}
        """
#
# index the alignments
# 
rule BamIndex:
    input:
        bam=rules.MapReads.output.align,
        #bam="alignments/align.{index}.bam",
    output:
        bai="alignments/align.{index}.bam.bai"
    params:
        cluster=" -pe serial 1 -l mfree=2G -l h_rt=24:00:00",
    shell:
        """
        samtools index {input.bam}
        """

#
# make a bed version of each bam file, this will be used to calculate coverage
#
rule BamToBed:
    input:
        bam=rules.MapReads.output.align,
        #bam="alignments/align.{index}.bam",
        bai=rules.BamIndex.output.bai,
        #bai="alignments/align.{index}.bam.bai"
    output:
        bed="alignments/align.{index}.bam.bed"
    params:
        cluster=" -pe serial 1 -l mfree=4G -l h_rt=24:00:00",
        samutils=config["mcutils"]
    shell:
        """
        samtools view {input.bam} | {params.samutils}/samToBed /dev/stdin --reportIdentity > {output.bed}
        """

#
# this rule creats a bed file that is incremented by 100 for every contig
# these will be the feautes upon which we calculate depth wtih bedtools
#
rule FaiToBed:
	input:
		asmfai=reference + ".fai",
	output:
		regions="coverage/regions.100.bed",
		regions1k="coverage/regions.1000.bed",
	params:
		cluster=" -pe serial 1 -l mfree=2G -l h_rt=24:00:00",
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
# turn the bed files into a bed files that have the coverage of each 100bp segment of the genome
# 
rule BedToCoverage:
    input:
        regions="coverage/regions.100.bed",
        bed=rules.BamToBed.output.bed,
		#bed="alignments/align.{index}.bam.bed",
    output:
        coverage="coverage/coverage.{index}.bed",
    params:
        cluster=" -pe serial 1 -l mfree=4G -l h_rt=24:00:00",
    shell:
        """
		# get coverage and then sort by contig and then pos
        bedtools coverage -a {input.regions} -b {input.bed} -mean | \
				sort -k1,1 -k2,2n > {output.coverage}
        #bedtools merge -i test.sort.bed -c 4 -o sum
        """
#
# merge the coverage for each 100bp region in the genome 
#
rule MergeBed:
	input:
		coverage=dynamic(rules.BedToCoverage.output.coverage),
		#coverage=dynamic("coverage/coverage.{index}.bed"),
	output:
		bedsort="coverage/all.merged.bed",
	benchmark:
		"benchmarks/merge.bed.txt"
	params:
		cluster=" -pe serial 1 -l mfree=32G -l h_rt=24:00:00",
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
# calculate coverage for 1k
#
rule GetOneKregionCoverage:
	input:
		bed="coverage/all.merged.bed",
		regions1k="coverage/regions.1000.bed",
	output:
		bed="coverage/all.1000.merged.bed",
	params:
		cluster=" -pe serial 1 -l mfree=16G -l h_rt=24:00:00",
	run:
		print("starting the read of all.merged.bed")
		cols = ["contig", "start", "end", "coverage"]
		dtypes = {"contig":str, "start":int, "end":int, "cov":float}
		bed = pd.read_csv(input["bed"], sep = "\t", header=None, names = cols, dtype=dtypes, engine="c")
		print("done reading")
		
		out = []
		for contig, group in bed.groupby("contig", as_index=False):
			pre = 0
			for nxt in range(10,len(group), 10):
				start = group.iloc[pre]["start"]
				end = group.iloc[nxt-1]["end"]
				coverage = group.iloc[pre:nxt]["coverage"].mean()
				out.append("{}\t{}\t{}\t{}\n".format(contig, start, end, coverage))
				pre = nxt
			nxt = len(group)	
			start = group.iloc[pre]["start"]
			end = group.iloc[nxt-1]["end"]
			coverage = group.iloc[pre:nxt]["coverage"].mean()
			out.append("{}\t{}\t{}\t{}\n".format(contig, start, end, coverage))
			print(contig)
		out = "".join(out)
		open(output["bed"], "w+").write(out)


#
# calculate coverage in 5K windows, and step 1K at a time
#
rule FiveKWindowStepOneK:
	input:
		bed="coverage/all.1000.merged.bed",
	output:
		bed="coverage/all.5000.merged.bed",
	params:
		cluster=" -pe serial 1 -l mfree=16G -l h_rt=24:00:00",
	run:
		cols = ["contig", "start", "end", "coverage"]
		dtypes = {"contig":str, "start":int, "end":int, "cov":float}
		bed = pd.read_csv(input["bed"], sep = "\t", header=None, names = cols, dtype=dtypes, engine="c")
		out = []
		for contig, group in bed.groupby("contig", as_index=False):
			length = len(group)
			df = group.as_matrix(columns = ["start", "end", "coverage"])
			for idx in range(length):	
				start = int(df[idx,0])
				endidx = idx + 5 
				if(endidx >= length):
					endidx = length 	
				end = int(df[endidx-1][1])
				coverage = df[idx:endidx, 2].mean()
				out.append("{}\t{}\t{}\t{}\n".format(contig, start, end, coverage))
			print( out[-1] )

		out = "".join(out)
		open(output["bed"], "w+").write(out)


#
# calcualte the average coverages 
#
rule GetCoverageStats:
	input:
		one="coverage/all.merged.bed",
		oneK="coverage/all.1000.merged.bed",
		fiveK="coverage/all.5000.merged.bed",
	output:
		one="coverage/all.stats.txt",
		oneK="coverage/all.1000.stats.txt",
		fiveK="coverage/all.5000.stats.txt",
	params:
		cluster=" -pe serial 1 -l mfree=16G -l h_rt=24:00:00",
	run:
		for infile, outfile in zip(input, output):
			bed = pd.read_csv(infile, sep = "\t", header=None,
					names=['contig', 'start', 'end',"coverage"])
			stats = bed["coverage"].describe()
			out = "mean_coverage\tstd_coverage\n{}\t{}\n".format(stats["mean"], stats["std"])
			open(outfile,"w+").write(out)

#
# make configuration files that have the mincov, maxcov, and mintotal 
#
rule GenerateMinAndMax:
	input:
		#stats="coverage/all.5000.stats.txt"
		stats="coverage/all.stats.txt"
	output:
		minmax="MinMax.sh",
		json="MinMax.json"
	params:
		cluster=" -pe serial 1 -l mfree=1G -l h_rt=24:00:00",
	run:
		stats = pd.read_csv(input["stats"], header = 0, sep = "\t")
		# the plust one is to make it round up
		maxcov = int(stats.iloc[0]["mean_coverage"] )
		mintotal=int(2.0*stats.iloc[0]["mean_coverage"] - 0.1*stats.iloc[0]["std_coverage"]+1)
		mincov = int(stats.iloc[0]["mean_coverage"] * 0.5 - 0.1*stats.iloc[0]["std_coverage"]+1)
		out = "export MINCOV={}\nexport MAXCOV={}\nexport MINTOTAL={}\n".format(mincov, maxcov, mintotal)
		open(output["minmax"], "w+").write(out)
		out2 = '{{\n\t"MINCOV" : {},\n\t"MAXCOV" : {},\n\t"MINTOTAL" : {}\n}}\n'.format(mincov, maxcov, mintotal)
		open(output["json"], "w+").write(out2)

#
# count the number of ovlapping bases between the repeate masking and all,merged,.bed
#
rule CountOverlappingRepeatElements:
	input:
		combined="coverage/all.merged.bed",
		#combined="coverage/all.5000.merged.bed",
		allR = "reference/all.repeats.bed",
	output:
		repeatCounted="coverage/all.repeatCounted.bed",
	params:
		cluster=" -pe serial 1 -l mfree=16G -l h_rt=24:00:00",
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
	if(row["repeatCount"] == 0):
		val = "+"
	return(val)

#
# confine bed file to only high coverage regions (hrc) and then merge any that are adj to eachother
#
rule BedForCollapses:
	input:
		combined="coverage/all.repeatCounted.bed",
		json="MinMax.json",
		fai=reference + ".fai",
	output:
		temp = "coverage/unmerged.100.collapses.bed",
		collapses="coverage/unmerged.collapses.bed",
	params:
		cluster=" -pe serial 1 -l mfree=16G -l h_rt=24:00:00",
	run:
		#cov = json.load(open(input["json"]))
		bed = pd.read_csv(input["combined"], sep = "\t", header=None, 
				names=['contig', 'start', 'end', 'coverage', "repeatCount"])

		# require high enough coverage
		stats = bed["coverage"].describe()
		mincov= 2.0*stats["mean"] - 0.1*stats["std"]
		print(mincov)
		bed = bed.ix[ bed["coverage"] >= mincov ]
		
		# marks the region as having or not having repeat content by strand
		bed["isNotRepeat"] = bed.apply(NoRepeatContent, axis=1)
		# writes to file before merging
		bed.to_csv(output["temp"], header=False, index=False, sep="\t" )
		# create a new file that merges adj regions that have the same strand (i.e. repeat content status) 
		shell("bedtools merge -i {output.temp} -d 2 -s -c 4,5 -o mean,mean > {output.collapses}")


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
		json="MinMax.json",
		fai=reference + ".fai",
		#collapses="marksHCR.bed",
		collapses="coverage/unmerged.collapses.bed",
		allR = "reference/all.repeats.bed",
	output:
		unf = "unfiltered.collapses.bed",
	params:
		cluster=" -pe serial 1 -l mfree=8G -l h_rt=24:00:00",
	run:	
		# read in the merged set
		HCR = pd.read_csv(input["collapses"], sep = "\t", header=None, 
				names=['contig', 'start', 'end', "notRepeat", 'coverage', "repeatPer"])
		# calcualte collapse length, +1 is because they are inclusive ranges on both sides
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

rule FilterCollapses:
	input:
		unf = "unfiltered.collapses.bed",
	output:
		collapses="collapses.bed",
		png ="SizeRepeatFilter.png",
	params:
		cluster=" -pe serial 1 -l mfree=8G -l h_rt=24:00:00",
	run:
		collapses = pd.read_csv(input["unf"], sep = "\t", header=None, 
				names=['contig', 'start', 'end', "notR", 'coverage', "RC", "length", "contigl", "distToEnd"])
		minsize = 9000
		maxRC = 75
		# plot what filter will be 
		cmd = "/net/eichler/vol2/home/mvollger/projects/abp/PlotFilterBySizeAndRepeatContent.R --bed {} --png {} --size {} --repeatContent {}".format(input["unf"], output["png"], minsize, maxRC)
		shell(cmd)
		# apply filter
		collapses = collapses.ix[(collapses["length"] >= minsize) & (collapses["RC"]<=maxRC)]
		
		outf = output["collapses"] 
		#outf = "temp.bed"
		collapses.to_csv(outf, header=False, index=False, sep="\t" )



# creates a regions file that has all the regions that are collapsed
# and creates a directory for eahc one of those regions
#
rule LocalAssembliesRegions:
	input:
		png ="SizeRepeatFilter.png",
		collapses="collapses.bed",
	output:
		regions="LocalAssemblies/regions.txt",
	params:
		cluster=" -pe serial 1 -l mfree=8G -l h_rt=24:00:00",
	run:
		df = pd.read_csv(input["collapses"], sep="\t", header=None)
		df["start"] = df[1]#/100
		df["start"] = df.start.map(int)
		df["end"] = df[2]#/100
		df["end"] = df.end.map(int)
		df["ID"] = df[0] + "." + df.start.map(str) + "." + df.end.map(str) + "/"
		df[["ID"]].to_csv(output["regions"], header=False, index=False, sep="\t")
		
		# for some reason making directories in a dynamic rule messes things up,
		# so i am going to make the collapse directories here
		rfile = open(output["regions"])
		dirsForShell = ""
		for line in rfile:
			region = "LocalAssemblies/" + line.strip()
			dirsForShell += region + " "
		rfile.close()
		# remove any old directories
		# shell('rm -rf LocalAssemblies/*.*.*')
		# add new direcotires 
		shell("mkdir -p " + dirsForShell)


#
# add a bed file to each region specifying where in asm they were 
#
rule LocalAssembliesBed:
	input:
		collapses="collapses.bed",
		regions="LocalAssemblies/regions.txt",
	output:
		bed=dynamic("LocalAssemblies/{region}/orig.bed"),
		rgn=dynamic("LocalAssemblies/{region}/orig.rgn"),
	params:
		cluster=" -pe serial 1 -l mfree=1G -l h_rt=1:00:00",
	run:
		rfile = open(input["regions"])
		regions=[]
		for line in rfile:
			region = "LocalAssemblies/" + line.strip()
			regions.append(region)
		rfile.close()

		# create reference and bed file
		cfile = open(input["collapses"])
		collapses = cfile.readlines()
		cfile.close()
		cmd = ""
		for line, region in zip(collapses,regions):
			line = line.split("\t")
			# for making the bed file
			bed = "{}\t{}\t{}".format(line[0], int(float(line[1])), int(float(line[2])) )
			cmd += "echo {} > {}/orig.bed; ".format(bed, region)
			# for making the regions file
			rgn = "{}:{}-{}".format(line[0], int(float(line[1])), int(float(line[2])) )
			cmd += "echo {} > {}/orig.rgn; ".format(rgn, region)

		shell(cmd)
#
# same as above but in a different format
# this rule requires a long latency wait time for some reason ~60 seconds 
#
'''
rule LocalAssembliesRgn:
	input:
		#bed=rule.LocalAssembliesBed.output.bed, #this does not work because the output is dynamic()
		bed = "LocalAssemblies/{region}/orig.bed"
	output:
		rgn = "LocalAssemblies/{region}/orig.rgn"
	params:
		cluster=" -pe serial 1 -l mfree=1G -l h_rt=1:00:00",
	run:
		myfile = open(input["bed"])
		bed = myfile.read().strip().split()
		myfile.close()
		assert len(bed) == 3
		rgn = "{}:{}-{}\n".format(bed[0], bed[1], bed[2])
		f = open(output["rgn"],"w+")
		f.write(rgn)
		f.close()
'''

#
# using the .rgn files make a fasta file consisting of the collapse 
#
rule LocalAssembliesRef:
	input:
		rgn="LocalAssemblies/{region}/orig.rgn",
		bed="LocalAssemblies/{region}/orig.bed",
		asm=config["asm"]
		#rgn=rules.LocalAssembliesRgn.output.rgn,
	output:
		refs="LocalAssemblies/{region}/ref.fasta",
	params:
		cluster=" -pe serial 1 -l mfree=4G -l h_rt=1:00:00",
	shell:
		"""
		region=$(cat {input.rgn})
		samtools faidx {input.asm} $region > {output.refs}
		"""
	'''
	run:
		region = open(input["rgn"]).read().strip()
		cmd = "samtools faidx " + input["asm"] + " " + region + " > " + output["refs"]
		shell(cmd)
	'''
#
# find the reads in all the bam files that map to that region
#
rule LocalAssembliesBam:
	input:
		refs=rules.LocalAssembliesRef.output.refs,
		rgn= rules.LocalAssembliesRef.input.rgn,
		bed= rules.LocalAssembliesRef.input.bed,
		#refs="LocalAssemblies/{region}/ref.fasta",
		#rgn="LocalAssemblies/{region}/orig.rgn",
	output:
		bams="LocalAssemblies/{region}/reads.orig.bam"
	params:
		cluster=" -pe serial 1 -l mfree=4G -l h_rt=24:00:00",
	run:
		import pysam
		myfile = open(input["rgn"])
		region = myfile.read().strip()
		myfile.close()

		bed = open(input["bed"]).read().strip()
		token = bed.split()
		
		allreads = None
		for idx, bam in enumerate(sorted(glob.glob("alignments/align.*.bam"))):
			samfile = pysam.AlignmentFile(bam, "rb")
			if(idx == 0 ):
				allreads = pysam.AlignmentFile(output["bams"], "wb", template=samfile)

			for read in samfile.fetch(token[0], int(token[1]), int(token[2])):
				allreads.write(read)
			samfile.close()		

			#out = output["bams"] + "." + str(idx) + ".bam"
			# not sure if it gets that only overlap
			#cmd = "samtools view -b {} {} > {}".format(bam, region, out)
			# failing becasue of a weird bin error with bedtools
			#cmd = "bedtools intersect -wa -abam {} -b {} > {}".format(bam, input["bed"], out)
			#shell(cmd)	
		allreads.close()
		
		#mergeCmd = "samtools merge {} {}.*.bam".format(output["bams"],  output["bams"])
		#shell(mergeCmd)
		# remove the temp files 
		#shell("rm LocalAssemblies/"+wildcards["region"]+"/reads.orig.bam.*.bam")

#
# copy ofver the min max stats so that ABP knows mincov, maxcov, mintotal
#
rule LocalAssembliesConfig:
	input:
		MinMax="MinMax.json",
		bams=rules.LocalAssembliesBam.output.bams,
		#bams="LocalAssemblies/{region}/reads.orig.bam"
	output:
		cov="LocalAssemblies/{region}/coverage.json",
	params:
		cluster=" -pe serial 1 -l mfree=1G -l h_rt=1:00:00",
		project=config["project"],
	run:
		f = open(input["MinMax"])
		out = ""
		for idx, line in enumerate(f.readlines()):
			if(idx == 1):
				out += '\t"project":"{}",\n'.format(params["project"]) 
			out += line
		open(output["cov"], "w+").write(out)
		'''
		cp {input.MinMax} {output.cov}
		sed -i '1 s/$/\n\t"project":"{params.project}",/' {output.cov}
		'''

#
# create the sequences in the reference that best match what I am generating 
#
GRCh38 = config["reference"]
fai = GRCh38 + ".fai"  
sa = GRCh38 + ".sa"  
#
# combine all of the ref.fastas so I can blasr them with one command, (much faster)
#
rule combineRefFasta:
	input:
		cov=dynamic(rules.LocalAssembliesConfig.output.cov),
		#cov=dynamic("LocalAssemblies/{region}/coverage.json"),
		dupref=dynamic(rules.LocalAssembliesRef.output.refs),
		#dupref=dynamic("LocalAssemblies/{region}/ref.fasta"),
	output:
		allref="LocalAssemblies/all.ref.fasta",
	params:
		cluster=" -pe serial 4 -l mfree=4G -l h_rt=3:00:00",
	shell:
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
rule duplicationsFasta1:
	input:
		dupref="LocalAssemblies/all.ref.fasta",
	output:
		dupsam="LocalAssemblies/all.ref.fasta.sam",
	params:
		blasr=config["blasr"],
		cluster=" -pe serial 4 -l mfree=16G -l h_rt=12:00:00",
	shell:
		"""
		# for some reason if i use the version of blasr in params.blasr it fails.
		# this is a bug I should probably share with mark
		# it looks like lots of alignments have 0 matching bases
		# so I will use the regualrly loaded blasr
		# {params.blasr}/alignment/bin/blasr  
		echo "here"	
		#blasr -clipping subread {input.dupref} {GRCh38} -nproc {threads} -sa {sa} \
		#		-sam -bestn 30 -minMatch 15 -maxMatch 20 \
		#		-out {output.dupsam}
		
		blasr -nproc 4 -sa {sa} -sam -out /dev/stdout \
			-minMatch 11 -maxMatch 20 -nCandidates 50 -bestn 30 \
			{input.dupref} {GRCh38} | \
			samtools view -h -F 4 - | samtools sort -m 4G -T tmp -o {output.dupsam}
		echo "here2"

		"""
		
#
# filter the sam file to only include high identity long contigs 
#
rule getHighIdentity:
	input:
		dupsam="LocalAssemblies/all.ref.fasta.sam",
	output:
		duptsv="LocalAssemblies/all.ref.fasta.identity.tsv",
	params:
		cluster=" -pe serial 1 -l mfree=1G -l h_rt=1:00:00",
	shell:
		"""
		#~mchaisso/projects/mcutils/bin/samToBed {input.dupsam} --reportIdentity | awk '{{ if ($3-$2 >8999 && $9 > 0.80) print $1":"$2"-"$3 }}' > {output.duptsv}
		~mchaisso/projects/mcutils/bin/samToBed {input.dupsam} --reportIdentity > {output.duptsv}
		"""


#
# generate two bed file for each assmebliy, one with the region of the collapse, and one with 100000 bp of slop 
# on either side
#
rule ConvertTsvToBedAndRgn:
	input:
		duptsv="LocalAssemblies/all.ref.fasta.identity.tsv",
	output:
		bedDone="LocalAssemblies/bed.done.txt",
		#dupbed=dynamic("LocalAssemblies/{region}/ref.fasta.long.bed")
	params:
		cluster=" -pe serial 1 -l mfree=1G -l h_rt=3:00:00",
	run:
		df = pd.read_csv( input["duptsv"], sep="\t", header=None,
				names=["contig", "start", "end", "read", "x", "y", "z", "z", "perID", "m", "mm", "i", "d"])
		df["length"] = df["end"] - df["start"]
		df=df.ix[ (df["length"]>=5000) & (df["perID"] > 0.80) ]
		df.reset_index(drop=True, inplace=True)
		
		allbed = "LocalAssemblies/all.ref.fasta.bed"
		df[["contig", "start", "end"]].to_csv(allbed, sep="\t", index=False, header=False)
		
		shell("bedtools slop -i {} -g {} -b 100000 > {}".format(allbed, fai, allbed+".slop"))	
		slop = pd.read_csv( allbed + ".slop", sep="\t", header=None, names=["contig", "start", "end"])
		df["longstart"] = slop["start"]
		df["longend"] = slop["end"]

		grouped = df.groupby(["read"])
		for name, group in grouped:
			match =  re.search('(.+):(\d+)-(\d+)/.*', name)
			if(match):
				region = "{}.{}.{}".format(match.group(1),match.group(2),match.group(3))
				outfile = "LocalAssemblies/{}/ref.fasta.bed".format(region)
				group[["contig", "start", "end"]].to_csv(outfile, sep="\t", index=False, header=False)
				outfile2 = "LocalAssemblies/{}/ref.fasta.long.bed".format(region)
				group[["contig", "longstart", "longend"]].to_csv(outfile2, sep="\t", index=False, header=False)
		shell("touch " + output["bedDone"])
					
#
# actaully fetch that region from the genome 
#
rule getReferenceSequences:
	input:
		#duplong="LocalAssemblies/{region}/ref.fasta.long.bed",
		bedDone="LocalAssemblies/bed.done.txt",
		regions="LocalAssemblies/regions.txt",
	output:
		#dup="LocalAssemblies/{region}/duplications.fasta",
		refDone="LocalAssemblies/README.txt",
	params:
		cluster=" -pe serial 1 -l mfree=4G -l h_rt=1:00:00",
	run:
		regions = open(input["regions"])
		for region in regions.readlines(): 
			region = region.strip()[:-1]
			bedfile = "LocalAssemblies/{}/ref.fasta.long.bed".format(region)
			fastafile = "LocalAssemblies/{}/duplications.fasta".format(region)
			cmd = "bedtools getfasta -fi {} -bed {} > {}".format(GRCh38, bedfile, fastafile)
			if(os.path.exists(bedfile)):
				shell(cmd)
		shell("touch " + output["refDone"])







