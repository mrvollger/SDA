import os
import numpy as np
import pandas as pd
import json
import re
import sys

#
# setup the env for each exacution 
#
SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)
snake_dir = SNAKEMAKE_DIR + "/"
shell.executable("/bin/bash")
shell.prefix("source %s/env_python3.cfg; " % SNAKEMAKE_DIR)


baxsets = {} 
f = open("input.fofn")
for line in f.readlines():
	line = line.strip()
	match = re.match(".*/(.+)\.[123]\.bax\.h5", line) 
	if(match is None):
		print("No match for {}".format(line), file = sys.stderr)
	else:
		batchid = match.group(1)
		#print(batchid)
		if(batchid not in baxsets):
			baxsets[batchid] = []
		baxsets[batchid].append(line)

# rmeovve items not in triplicate 
for key in baxsets:
	values = baxsets[key]
	if(len(values) != 3):
		print("Not using:{}".format(values), file = sys.stderr)
		del baxsets[key]

baxids = baxsets.keys()
#print(baxsets)

# get baxes involved with a set
def getBaxs(baxid):
	return(baxsets[str(baxid)])


localrules: all, final, fofn,


rule all:
	input:
		"final"


rule Bax2Bam:
	input:
		baxset = getBaxs,
	output:
		bam = "bams/{ID}.subreads.bam",
		pbi = "bams/{ID}.subreads.bam.pbi",
		scrap = temp("bams/{ID}.scraps.bam"),
		scrapindex = temp("bams/{ID}.scraps.bam.pbi"),
	params:
		cluster = " -l mfree=8G "
	run:
		prefix = output["bam"][:-13]
		shell("bax2bam {} -o {}".format(input["baxset"], prefix) )
		


rule ProtectAndMd5:
	input:
		bam = "bams/{ID}.subreads.bam",
		pbi = "bams/{ID}.subreads.bam.pbi",
	output:
		md5 = temp("bams/{ID}.md5"),
	params:
		cluster = " -l mfree=8G "
	run:
		# write protect file
		shell("chmod g-w {input.bam}")
		shell("chmod g-w {input.pbi}")

		# add to md5sum 
		shell("md5sum {input.bam} {input.pbi} > {output.md5}")


rule fofn:
	input:	
		bam = expand("bams/{ID}.subreads.bam", ID=baxids),
		md5 = expand("bams/{ID}.md5", ID=baxids),
	output:
		fofn = protected("bams.fofn"),
		md5 = protected("md5.txt"),
	run:
		out = ""
		for bam in input["bam"]:
			out += os.path.abspath(bam) + "\n"
		open(output["fofn"], "w+").write(out)
		
		# make merged md5
		shell("cat {input.md5} > {output.md5}")


rule final:
	input:
		fofn = "bams.fofn",
	output:
		final = temp("final")
	shell:
		"touch {output.final}"


