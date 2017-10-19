import os
import tempfile
import numpy as np
import pandas as pd
import json
import glob
from pprint import pprint

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



localrules: all, MergeCoverage, FaiToBed, getCoverageStats, 
			GenerateMinAndMax, bedForCollapses, MergeBed, 
			LocalAssembliesBed, LocalAssembliesRegions,
			LocalAssembliesRef,
			LocalAssembliesBam,


rule all:
	input:
		combined="coverage/all.sorted.merged.bed",
		stats="coverage/all.stats.txt",
		minmax="MinMax.sh",
		collapses="collapses.bed",
		regions="LocalAssemblies/regions.txt",
		bams=dynamic("LocalAssemblies/{region}/reads.orig.bam"),
		#asmSA=config["asm"] + ".sa",
		#regions="coverage/regions.bed",
		#split=dynamic("reads.{index}.fofn"),
		#align=dynamic("alignments/align.{index}.bam"),
		#alignbed=dynamic("alignments/align.{index}.bam.bed"),
		#coverage=dynamic("coverage/coverage.{index}.bed"),

#
# If a suffix arry does not yet exist for the assembly, build this.
#
rule MakeASM_SA:
    input:
        asm=config["asm"]
    output:
        asmsa=config["asm"] + ".sa"
    params:
        sge_opts=" -pe serial 1 -l mfree=48G -l h_rt=6:00:00",
        blasr=config["blasr"]
    shell:
        """
        {params.blasr}/alignment/bin/sawriter {input.asm}
        """

rule SplitFOFN:
    input:
        fofn=config["reads"]
    output:
        fofnSplit=dynamic("reads.{index}.fofn")
    params:
        sge_opts=" -pe serial 1 -l mfree=1G -l h_rt=1:00:00",
        baxPerJob=config["bax_per_job"]
    shell:
        """
        split --numeric --lines {params.baxPerJob} {input.fofn} reads.
        for f in `ls reads*`; do
        mv $f $f.fofn
        done
        """

#
#  For read depth, and other future steps, it is necessary to map reads back to the assembly.
#
rule MapReads:
    input:
        asm=config["asm"],
        asmsa=config["asm"] + ".sa",
        fofn="reads.{index}.fofn"
    output:
        align="alignments/align.{index}.bam"
    params:
        sge_opts=" -pe serial 8 -l mfree=8G -l h_rt=24:00:00",
        blasr=config["blasr"]
    shell:
        """
        mkdir -p alignments
        {params.blasr}/alignment/bin/blasr {input.fofn} {input.asm} \
            -nproc 8 -out /dev/stdout -minMapQV 0 -bestn 2 -advanceExactMatches 15 -sam | \
            samtools view -bS - | \
            samtools sort -T $TMPDIR/blasr -@ 8 -m4G - -o {output.align}
        samtools index {output.align}
        """
        

rule BamToBed:
    input:
        bam="alignments/align.{index}.bam"
    output:
        bed="alignments/align.{index}.bam.bed"
    params:
        sge_opts=" -pe serial 1 -l mfree=2G -l h_rt=24:00:00",
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
        asmfai=config["asmfai"],
    output:
        regions="coverage/regions.bed",
    #params:
        #sge_opts=" -pe serial 1 -l mfree=2G -l h_rt=24:00:00",
    run:
        shell("mkdir -p coverage")
        fai = open(input["asmfai"])
        out = ""
        for line in fai:
            token = line.strip().split("\t")
            length = int(token[1])
            contig = token[0]
            cur = 0
            for nxt in range(100, length, 100):
                out += "{}\t{}\t{}\n".format(contig, cur, nxt-1)
                cur = nxt
            out += "{}\t{}\t{}\n".format(contig, cur, length-1)
            #print(out)
            #print("{}\t{}".format(token[0], token[1]))
        outfile = open(output["regions"], "w+")
        outfile.write(out)



rule BedToCoverage:
    input:
        regions="coverage/regions.bed",
        bed="alignments/align.{index}.bam.bed",
    output:
        coverage="coverage/coverage.{index}.bed",
    params:
        sge_opts=" -pe serial 1 -l mfree=4G -l h_rt=24:00:00",
    shell:
        """
		# get coverage and then sort by contig and then pos
        bedtools coverage -a {input.regions} -b {input.bed} -mean | \
				sort -k1,1 -k2,2n > {output.coverage}
        #bedtools merge -i test.sort.bed -c 4 -o sum
        """
#
# merge the coverage files using bedtools
#
rule MergeBed:
	input:
		coverage=dynamic("coverage/coverage.{index}.bed"),
	output:
		sorted="coverage/all.sorted.bed",
	shell:
		"""
		#  contig, pos, merge the input files
		echo $TMPDIR
		sort -T $TMPDIR -k1,1 -k2,2n -S 8G -m {input.coverage} > {output.sorted}
		"""

#
# merge the coverage files using bedtools
#
'''
rule PythonMergeBed:
	input:
		coverage=dynamic("coverage/coverage.{index}.bed"),
	output:
		sorted="coverage/all.sorted.bed",
	run:
'''		

rule MergeCoverage:
	input:
		sorted="coverage/all.sorted.bed",
	output:
		combined="coverage/all.sorted.merged.bed",
	shell:
		"""
		# our output should already be sorted so we can use merge straight off
		bedtools merge -i {input.sorted} -c 4 -o sum > {output.combined}
		"""

#
# calcualte the average coverages 
#
rule getCoverageStats:
	input:
		combined="coverage/all.sorted.merged.bed",
	output:
		stats="coverage/all.stats.txt"
	run:
		bed = open(input["combined"])
		covs = []
		for line in bed:
			token = line.strip().split("\t")
			covs.append(float(token[3]))
		covs = np.array(covs)
		mean = np.mean(covs)
		median = np.median(covs)
		std = np.std(covs)
		out = "mean_coverage\tmedian_coverage\tstd_coverage\n{}\t{}\t{}\n".format(mean, median, std)
		print(out)
		open(output["stats"],"w+").write(out)



rule GenerateMinAndMax:
	input:
		stats="coverage/all.stats.txt"
	output:
		minmax="MinMax.sh",
		json="MinMax.json"
	run:
		stats = pd.read_csv(input["stats"], header = 0, sep = "\t")
		maxcov = int(stats.iloc[0]["median_coverage"])
		mintotal=int(maxcov*1.5)
		mincov = int(maxcov * 0.6 )
		out = "export MINCOV={}\nexport MAXCOV={}\nexport MINTOTAL={}\n".format(mincov, maxcov, mintotal)
		open(output["minmax"], "w+").write(out)
		out2 = '{{\n\t"MINCOV" : {},\n\t"MAXCOV" : {},\n\t"MINTOTAL" : {}\n}}\n'.format(mincov, maxcov, mintotal)
		open(output["json"], "w+").write(out2)


rule bedForCollapses:
	input:
		combined="coverage/all.sorted.merged.bed",
		json="MinMax.json",
		fai=config["asmfai"],
	output:
		temp = "coverage/unmerged.collapses.bed",
		collapses="collapses.bed",
	run:
		cov = json.load(open(input["json"]))
		bed = pd.read_csv(input["combined"], sep = "\t", header=None, names=['contig', 'start', 'end', 'coverage'])
		# get coverage stats per contig, this will be used later for filtering 
		coverageByContig = bed.groupby('contig')['coverage'].describe()
	 	# the index is the contig, going to make it a colunm so I can merge
		coverageByContig['contig'] = coverageByContig.index
		
		# require high enough coverage
		bed = bed.ix[ bed["coverage"] >= cov["MINTOTAL"] ]
		bed.to_csv(output["temp"], header=False, index=False, sep="\t" )
		# create a new file that has regions within 1001 base pairs of one another merged, changed -d 1001 to -d 1
		shell("bedtools merge -i {output.temp} -d 1 -c 4 -o mean > {output.collapses}")
	
		# read in the merged set
		merged = pd.read_csv(output["collapses"], sep = "\t", header=None, names=['contig', 'start', 'end', 'coverage'])
		# read in the length of the contigs
		fai = pd.read_csv(input["fai"], sep = "\t", header=None, names=['contig', 'length', 'OFFSET', 'LINEBASES', "LINEWIDTH"])
		fai = fai[["contig", "length"]]
		
		# this adds the contigs length to the collapse 
		merged = pd.merge(merged, fai, on='contig', how='outer')
		# require the collapse to be within 50000 bp of the ned of a contig
		merged["distFromEnd"]=pd.concat([merged["start"], merged["length"]-merged["end"]], axis=1).min(axis=1)
		merged = merged.ix[merged["distFromEnd"] <= 50000]
		# require the contig be over a certain length, 100000
		merged = merged.ix[ merged["length"] >= 300000 ]


		# merge the previously calcualted contig stats into each row 
		merged = pd.merge(merged, coverageByContig, on='contig', how ='outer' )

		# calcualte collapse length	
		merged["clength"] = merged["end"] - merged["start"]
		merged=merged.ix[merged["coverage"]>=cov["MINTOTAL"]]
		# require a collpase over 10000 base pairs
		merged=merged.ix[merged["clength"]>=10000]
		# write to file
		merged.to_csv(output["collapses"], header=False, index=False, sep="\t" )



rule LocalAssembliesRegions:
	input:
		collapses="collapses.bed",
	output:
		#refs=dynamic("LocalAssemblies/{region}/ref.fasta"),
		regions="LocalAssemblies/regions.txt",
	run:
		df = pd.read_csv(input["collapses"], sep="\t", header=None)
		df["start"] = df[1]/100
		df["start"] = df.start.map(int)
		df["end"] = df[2]/100
		df["end"] = df.end.map(int)
		df["ID"] = df[0] + "_" + df.start.map(str) + "_" + df.end.map(str) + "/"
		df[["ID"]].to_csv(output["regions"], header=False, index=False, sep="\t")



rule LocalAssembliesBed:
	input:
		collapses="collapses.bed",
		regions="LocalAssemblies/regions.txt",
	output:
		bed=dynamic("LocalAssemblies/{region}/orig.bed"),
		rgn=dynamic("LocalAssemblies/{region}/orig.rgn"),
	run:
		rfile = open(input["regions"])
		dirsForShell = ""
		regions=[]
		for line in rfile:
			region = "LocalAssemblies/" + line.strip()
			regions.append(region)
			dirsForShell += region + " "
		shell("mkdir -p " + dirsForShell)
		
		# create reference and bed file
		cfile = open(input["collapses"]).readlines()
		for line, region in zip(cfile,regions):
			line = line.split("\t")
			seq = "{}:{}-{}".format(line[0], int(float(line[1])), int(float(line[2])) )
			cmd = "echo " + "{}\t{}\t{}".format(line[0], int(float(line[1])), int(float(line[2])))+" > "+region+"orig.bed;"
			cmd += "echo " + seq + " > " + region + "orig.rgn;"
			shell(cmd)


rule LocalAssembliesRef:
	input:
		rgn="LocalAssemblies/{region}/orig.rgn",
		asm=config["asm"],
	output:
		refs="LocalAssemblies/{region}/ref.fasta",
	#params:
	#	sge_opts=" -pe serial 1 -l mfree=4G -l h_rt=1:00:00",
	run:
		region = open(input["rgn"]).read().strip()
		cmd = "samtools faidx " + input["asm"] + " " + region + " > " + output["refs"]
		shell(cmd)

rule LocalAssembliesBam:
	input:
		refs="LocalAssemblies/{region}/ref.fasta",
		rgn="LocalAssemblies/{region}/orig.rgn",
	output:
		bams="LocalAssemblies/{region}/reads.orig.bam"
	#params:
	#	sge_opts=" -pe serial 1 -l mfree=4G -l h_rt=24:00:00",
	run:
		region = open(input["rgn"]).read().strip() 
		for idx, bam in enumerate(glob.glob("alignments/align.*.bam")):
			out = output["bams"] + "." + str(idx) + ".bam"
			cmd = "samtools view -b {} {} > {}".format(bam, region, out)
			shell(cmd)	

		mergeCmd = "samtools merge {} {}.*.bam".format(output["bams"],  output["bams"])
		shell(mergeCmd)
		# remove the temp files 
		shell("rm LocalAssemblies/*/reads.orig.bam.*.bam")



