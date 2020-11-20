import os
import sys
import tempfile 
import numpy as np
import re
import pysam
from collections import Counter
from Bio import SeqIO
import pandas as pd
import tempfile 

snake_dir = os.path.dirname(workflow.snakefile) + "/"
CWD = os.getcwd()

python3 =  f"cd {snake_dir} && source env_sda.sh && cd {snake_dir}envs/ && source env_python3.sh && cd {CWD}" 

shell.executable("/bin/bash")
shell.prefix(f"{python3} && set -eo pipefail; ")

base = snake_dir + "scripts/"

TRF = snake_dir + "/externalRepos/trf-v4.09/bin/trf"

#
# Get tmp dir
#
if "TMPDIR" in os.environ:
    TMPDIR = os.environ['TMPDIR']
elif "TMPDIR" in config:
    TMPDIR = config['TMPDIR']
else:
    TMPDIR = tempfile.gettempdir()


#
# script locations and configurations 
#
INPUT = config["input"]
BAM=False
if(INPUT[-4:]==".bam"):
  BAM=INPUT
  assert os.path.exists(BAM+".bai"), "Input bam must be indexed!"
else:
  FOFN = INPUT
  READS =  [ line.strip() for line in open(FOFN).readlines() ] 
  IDS = list(range(len(READS ) ) )

PLAT = config["platform"].upper()
REF = os.path.abspath( config["ref"] )
MINALN = config["minaln"] 
BANDWIDTH = config["bandwidth"] 
PRE	= config["prefix"]
DIR	= config["dir"]
LRT = config["lrt"]
MINNUMSHARED = config["minNumShared"]
MAXPOSREP = config["maxPosRep"]
MINCUTSIZE = config["minCutSize"]
MINCUTLEN = config["minCutLen"]
ITERATIONS = config["iterations"]
RM_DB = config["species"]
ASSEMBLERS = config["assemblers"]
DEBUG = config["debug"]
MAX_TIME="300m"
# window size of calcualting coverage
WINDOW = 1000
# minimum size for a collapse
MINCOLLEN = 15000
# maximum percentage of a collapse that can be common repeat 
MAXCR = 75

# get the reference size
REF_GB=os.stat(REF).st_size/(1024*1024*1024)
I_G=4
if(REF_GB > 3.5):
  sys.stderr.write("You have a large reference! Increasing mem for alignment.")
  I_G = int(REF_GB) + 1


# make sure index is yonger
FAI = REF + ".fai"
if os.path.getmtime(FAI) < os.path.getmtime(REF):
  shell("samtools faidx {REF}")

#
# geting ready to run TRF and RM by splitting up the genome into multiple parts to run seperatly 
#
splitSize = 200
recs = open(FAI).readlines()
if(splitSize > len(recs)):
	splitSize = len(recs)
FRACS = list(range(splitSize))


#
# function to remove temp files
#
def tempd(File):
	if(DEBUG):
		return(File)
	return(temp(File))


wildcard_constraints:
	PRE = PRE,
	DIR = DIR,
	ID = "\d+",
	FRAC = "\d+",

onsuccess:
	sys.stderr.write("SDA DENOVO FINISHED\n")
onerror:
	sys.stderr.write("SDA DENOVO FAILED\n")

rule all:	
	input:
		hcr = f"{DIR}/coverage/{PRE}.collapses.bed",
		final = f"{DIR}/{PRE}.done",
		bam=f"{DIR}/{PRE}.reads.bam",
		bai=f"{DIR}/{PRE}.reads.bam.bai",


###################################################################################
#                                                                                 #
#             ALIGN READS TO THE INPUT REFERENCE                                  #
#                                                                                 #
###################################################################################

def get_reads(wildcards):
  ID = int(str(wildcards.ID))
  f = READS[ID]
  if PLAT in ["ONT"]:
    sys.stderr.write("Testing input type\n")
    assert re.match(".+\.(fa|fasta|fq|fastq)(\.gz)?", f), f +" does not match: .+\.(fa|fq|fasta|fastq)(\.gz)?"
  return(f)

if(PLAT in ["CCS", "SUBREAD"] ):
	rule pbmm2:
		input:
			reads = get_reads,
			ref = REF,
		output:
			tempd("{DIR}/{PRE}.{ID}.reads.bam")
		resources:
			mem=4+I_G
		threads: 8
		shell: """ 
# SUBREAD gives bettern alignmetns for CCS reads than CCS
pbmm2 align -j {threads} \
	--preset SUBREAD -N 50  --min-length {MINALN} -r {BANDWIDTH} \
	--sample FAKE_SAMPLE \
	{input.ref} {input.reads} | \
	samtools view -u -F 2308 - | \
	samtools sort -@ {threads} -m 4G -T {TMPDIR}/pbmm2 -o {output}
"""

elif( PLAT in ["ONT"] ): 
	rule minimap2:
		input:
			reads = get_reads,
			ref = REF,
		output:
			tempd("{DIR}/{PRE}.{ID}.reads.bam")
		resources:
			mem=4+I_G
		threads: 8
		shell:"""
minimap2 \
  -I {I_G}G \
	-ax map-ont \
	--eqx -L \
	-R '@RG\\tID:MINIMAP\\tSM:FAKE_SAMPLE\\tPL:ONT' \
	-t {threads} \
	-m {MINALN} -r {BANDWIDTH} \
	{input.ref} {input.reads} | \
	samtools view -u -F 2308 - | \
	samtools sort -@ {threads} -m 4G -T {TMPDIR}/pbmm2 -o {output}"""

else:
	sys.stderr.write("Platform {} not recongnized!\n".format(PLAT))

if(BAM):
  rule merge_bam:
    input:
      bam = BAM,
    output:
      bam = "{DIR}/{PRE}.reads.bam",
    resources:
      mem=1
    threads: 1
    shell:"""
  ln -s $(readlink -f {input.bam}) {output.bam}
  """	
  rule index_bam:
    input:
      bam=rules.merge_bam.output.bam,
      bai=BAM+".bai",
    output:
      bai="{DIR}/{PRE}.reads.bam.bai"
    resources:
      mem=1
    threads: 1
    shell:"""
  ln -s $(readlink -f {input.bai}) {output.bai}
  """
else:
  rule merge_bam:
    input:
      bams = expand("{{DIR}}/{{PRE}}.{ID}.reads.bam", ID=IDS),
    output:
      bam = "{DIR}/{PRE}.reads.bam",
    resources:
      mem=4
    threads: 12
    shell:"""
  samtools merge -@ {threads} {output.bam} {input.bams}
  """	

  rule index_bam:
    input:
      bam=rules.merge_bam.output.bam,
    output:
      bai="{DIR}/{PRE}.reads.bam.bai"
    resources:
      mem=16
    threads: 1
    shell:"""
  samtools index {input}
  """

#
# this rule creats a bed file that is incremented by 1000 for every contig
# these will be the feautes upon which we calculate depth wtih bedtools
#
rule fai_to_bed:
  input:	
    asmfai= REF + ".fai",
  output:
    regions=tempd("{DIR}/coverage/{PRE}.regions.bed"),
  resources:
    mem=16
  threads: 1
  shell:"""
bedtools makewindows -g {input.asmfai} -w {WINDOW} > {output.regions}
"""

rule bam_to_coverage:
	input:
		bam=rules.index_bam.input.bam,
		#bai=rules.index_bam.output.bai,
		genome=FAI,
		regions=rules.fai_to_bed.output.regions,
	output:
		cov=tempd("{DIR}/coverage/{PRE}.coverage.bed"),
	resources:
		mem=16
	threads: 1
	shell:"""
# get coverage and then sort by contig and then pos
bedtools coverage -bed -mean -sorted \
           -g {input.genome} -a {input.regions} \
           -b <( samtools view -b -F 2308 {input.bam} ) | \
	sort -k 1,1 -k2,2n > {output.cov}
"""




###################################################################################
#                                                                                 #
#             TRF AND REPEATMASKER ON THE INPUT REFERENCE                         #
#                                                                                 #
###################################################################################
MAX_RM_LEN = 48
rule split_ref:
	input:
		ref=REF,
	output:
		split = tempd(expand("{{DIR}}/common_repeats/{{PRE}}.ref.{FRAC}.fasta", FRAC=FRACS) ),
	resources:
		mem=16
	threads: 1
	run:
		seqs = list(SeqIO.parse(input["ref"], "fasta"))
		toWrite = {}
		count = 0
		for idx, seq in enumerate(seqs):
			if(count not in toWrite):
				toWrite[count] = []

			seq.id = str(idx)
			seq.name = str(idx)
			seq.description = str(idx)

			toWrite[count].append(seq)
			count += 1
			if(count == splitSize):
				count = 0

		for key in toWrite:
			print(key, len(toWrite[key]))
			SeqIO.write(toWrite[key], output["split"][key], "fasta")

#
# Run RepeatMasker
#
rule RepeatMasker:
  input:
    split = "{DIR}/common_repeats/{PRE}.ref.{FRAC}.fasta",
  output:
    out = tempd("{DIR}/common_repeats/{PRE}.ref.{FRAC}.fasta.out"),
    tbl = tempd("{DIR}/common_repeats/{PRE}.ref.{FRAC}.fasta.tbl"),
    #cat = tempd("{DIR}/common_repeats/{PRE}.ref.{FRAC}.fasta.cat"),
    #ref = tempd("{DIR}/common_repeats/{PRE}.ref.{FRAC}.fasta.ref"),
    #msk = tempd("{DIR}/common_repeats/{PRE}.ref.{FRAC}.fasta.masked"),
  log:
    "{DIR}/logs/RepeatMasker.{PRE}.ref.{FRAC}.stdout.log"
  resources:
    mem=8,
  threads:4
  run:
    rmdir = os.path.dirname(input["split"])
    rm_lib = f"-species {RM_DB}"
    if( os.path.exists(RM_DB)):
      rm_lib = f"-lib {RM_DB}"
    shell("""		
RepeatMasker \
    {rm_lib} \
    -e ncbi \
    -dir {rmdir} \
    -pa {threads} \
    {input.split} > {log}

# clean outputfiles that only sometimes exist depending on RM version
rm -f \
  {wildcards.DIR}/common_repeats/{wildcards.PRE}.ref.{wildcards.FRAC}.fasta.masked* \
  {wildcards.DIR}/common_repeats/{wildcards.PRE}.ref.{wildcards.FRAC}.fasta.cat* \
  {wildcards.DIR}/common_repeats/{wildcards.PRE}.ref.{wildcards.FRAC}.fasta.ref* 
    """)

#
# Run TRF
#
rule TRF:
	input:
		split = "{DIR}/common_repeats/{PRE}.ref.{FRAC}.fasta",
	output:
		trf = tempd("{DIR}/common_repeats/{PRE}.ref.{FRAC}.trf.dat"),
	resources:
		mem=16,
	threads:1
	run:
		pre = os.path.basename(input["split"])
		fasta = os.path.abspath(input["split"])
		dat = os.path.abspath(output["trf"])
		trfdir = os.path.dirname(input["split"])
		param = ["2", "7", "7", "80", "10", "50", "2000"]
		trfparam = " ".join(param)
		trfout = ".".join(param)
		# adding max runtime to trf because sometimes it stalls forever. 
		shell("""timeout {MAX_TIME} {TRF} {fasta} {trfparam} -h -ngs > {output.trf} || touch {output.trf} """ )



rule merge_trf_rm:
	input:
		rm = expand(f"{DIR}/common_repeats/{PRE}.ref.{{FRAC}}.fasta.out", FRAC = FRACS),
		asmfai= REF + ".fai",
		trf = expand(f"{DIR}/common_repeats/{PRE}.ref.{{FRAC}}.trf.dat", FRAC = FRACS),
	output:
		crtmp = temp(f"{DIR}/common_repeats/{PRE}.common_repeats.bed"),
		cr = f"{DIR}/common_repeats/{PRE}.common_repeats.sort.merge.bed",
		rm = f"{DIR}/common_repeats/{PRE}.rm.all.tbl",
		trf = f"{DIR}/common_repeats/{PRE}.trf.all.tbl",
	resources:
		mem=32,
	threads:1
	run:
		fai = pd.read_csv(input["asmfai"], sep="\t", names=["contig", "len", "x", "y", "z"] )
		convert = { idx:contig for idx, contig in enumerate(fai["contig"]) }
		

		#
		# PARSE TRF
		#
		trfnames = 'contig start end PeriodSize CopyNumber ConsensusSize PercentMatches PercentIndels Score A C G T Entropy Motif Sequence'.split()	
		trf= []
		for ftrf in input["trf"]:
			chrom = None
			sys.stderr.write(ftrf + "\n" )
			with open(ftrf, 'r') as dat:
				for line in dat:
					splitline = line.split()
					if( line.startswith("Sequence:") ):
						chrom = int(line.split()[1].strip())
						#sys.stderr.write(chrom + "\n")
					elif( line.startswith("@") ):
						chrom = int(splitline[0][1:].strip()) # grab everything after the @ in the first word
						#sys.stderr.write(chrom + "\n")
					else:
						# Catch index errors when line is blank
						try:
							# Check if in header sequence (all non-header lines start with an int: start pos)
							try:
								int(splitline[0])
							except ValueError:
								continue
							trf.append([chrom] + splitline[ 0: (len(trfnames)-1) ] )
						except IndexError:
							pass
		trf = pd.DataFrame(trf, columns=trfnames)
		print(trf.shape )
		
		trf["start"] = trf["start"].astype(int)
		trf.sort_values(by=["contig", "start"], inplace=True)
		print("done sorting trf")
		
		# insert original contig names
		#trf["contig"].replace(convert, inplace=True) # This function is very slow for no good reason.
		trf["contig"] = trf["contig"].map(convert.get)
		print("done converting trf")
		
		trf.to_csv(output["trf"], sep="\t", index=False)
		print("done writing trf")

		#
		# PARSE REPEAT MASKER
		#
		
		rms = []
		for frm in input["rm"]:
			sys.stderr.write(frm + "\n" )
			rms.append( pd.read_csv(frm, delim_whitespace=True, header=None, skiprows=[0,1,2], comment="*",
				names = ["score", "div", "del", "ins", "contig", "start", "end",
					 "q_left", "strand", "repeat", "class", "r_st", "r_en", "r_left", "id"]) )
		
		rm = pd.concat(rms, ignore_index=True)
		print(rm.shape )

		rm.sort_values(by=["contig", "start"], inplace=True)
		print("done sorting rm")

		# insert original contig names
		#rm["contig"].replace(convert, inplace=True)
		rm["contig"] = rm["contig"].map(convert.get)
		print("done converting rm")
		
		rm.to_csv(output["rm"], sep="\t", index=False)
		print("done writing rm")
	

		#
		# WRITE TO BED
		#
		bed = ["contig", "start", "end"]
		cm = pd.concat([rm[bed], trf[bed]], ignore_index=True)
		cm.to_csv(output["crtmp"], sep="\t", header=False, index=False)
		shell("bedtools sort -i {output.crtmp} | bedtools merge -i - > {output.cr}")



###################################################################################
#                                                                                 #
#             FIND HIGH COVERAGE REGIONS                                          #
#                                                                                 #
###################################################################################
rule simple_high_coverage_regions:
    input:
        cov = rules.bam_to_coverage.output.cov,
        asmfai= REF + ".fai",
    output:
        stats = "{DIR}/coverage/{PRE}.simple.coverage.stats",
        collapse = "{DIR}/coverage/{PRE}.simple.collapses.bed",
    resources:
        mem=16,
    threads:1
    run:
        bed = pd.read_csv( input["cov"], sep = "\t", header=None, names=['contig', 'start', 'end',"coverage"])
# I want to eliminte the really low or really high coverage things because they are probably
# not assembled correctly and then assecess what the mean and standard deviation is
        top = bed.coverage.quantile(.90)
        bot = bed.coverage.quantile(.10)

# save stats like mean coverage 
        stats = bed["coverage"][ (bed.coverage < top) & ( bed.coverage > bot) ].describe()
        out = "mean_coverage\tstd_coverage\n{}\t{}\n".format(stats["mean"], stats["std"])
        open(output["stats"], "w+").write(out)

# filter for high coverage regsion
        MINCOV = stats["mean"]  + 3 * np.sqrt(stats["mean"])
        shell("""
awk '{{ if ($4 > {MINCOV}) print;}}' {input.cov} \
    | bedtools merge -d 10 -c 4,4 -o mean,median \
    | awk '{{ if ($3-$2 > {MINCOLLEN}) {{ print $0"\t"$3-$2;}} }}' > {output.collapse} """)

		  

rule count_cm_per_window:
	input:
		cr = rules.merge_trf_rm.output.cr,
		cov = rules.bam_to_coverage.output.cov,
	output:
		cov = "{DIR}/coverage/{PRE}.coverage.repeat_counted.bed",
	resources:
		mem=16,
	threads:1
	shell:"""
# count number of overlaping bases with cm | eliminate extra colums | merge overlapping entries and calculate sum 
bedtools intersect -a {input.cov} -b {input.cr} -wao | cut -f 1,2,3,4,8 | bedtools merge -c 4,5 -o mean,sum > {output.cov} 
# final output:
# contig\tstart\tend\tcoverage\tcommon repeat bases
"""


rule high_coverage_regions:
	input:
		cov=rules.count_cm_per_window.output.cov, 
		asmfai= REF + ".fai",
	output:
		stats = "{DIR}/coverage/{PRE}.coverage.stats",
		hcr = "{DIR}/coverage/{PRE}.collapses.bed",
		hcr_cm = "{DIR}/coverage/{PRE}.collapses.with.cm.bed",
	resources:
		mem=16,
	threads:1
	run:
		bed = pd.read_csv( input["cov"], sep = "\t", header=None, names=['contig', 'start', 'end',"coverage", "cr"])

		# I want to eliminte the really low or really high coverage things because they are probably
		# not assembled correctly and then assecess what the mean and standard deviation is
		top = bed.coverage.quantile(.90)
		bot = bed.coverage.quantile(.10)
		
		# save stats like mean coverage 
		stats = bed["coverage"][ (bed.coverage < top) & ( bed.coverage > bot) ].describe()
		out = "mean_coverage\tstd_coverage\n{}\t{}\n".format(stats["mean"], stats["std"])
		open(output["stats"], "w+").write(out)
		

		# filter for high coverage regsion
		MINCOV = stats["mean"]  + 3 * np.sqrt(stats["mean"])
		shell("""
awk '{{ if ($4 > {MINCOV}) print;}}' {input.cov} \
	| bedtools merge -d 10 -c 4,4,5 -o mean,median,sum \
	| awk '{{ if ($3-$2 > {MINCOLLEN}) {{ print $0"\t"$3-$2;}} }}' > {output.hcr_cm} """)

		shell(""" 
awk '{{ if ($6/$7*100 <= {MAXCR}) {{ print $0}} }}' {output.hcr_cm} > {output.hcr} """)
		
###################################################################################
#                                                                                 #
#             RUN SDA ON HIGH COVERAGE REGIONS                                    #
#                                                                                 #
###################################################################################

COL_DIR = "{DIR}/{PRE}.LocalAssemblies/region_{LA_ID}"
COL_RGN_FMT = os.path.join( COL_DIR , "rgn.bed")
COL_REF_FMT = os.path.join( COL_DIR , "ref.fasta")
COL_BAM_FMT = os.path.join( COL_DIR , "reads.orig.bam")
COL_SDA_FMT = os.path.join( COL_DIR , "region_{LA_ID}.done")

checkpoint local_asm_dirs:
	input:
		ref = REF,
		hcr = rules.high_coverage_regions.output.hcr,
	output:
		LAs = directory("{DIR}/{PRE}.LocalAssemblies/")
	resources:
		mem=16,
	threads:1
	run:
		for LA_ID, line in enumerate(open(input["hcr"]).readlines()):
			rgn = COL_RGN_FMT.format(DIR=DIR, PRE=PRE, LA_ID=LA_ID)
			shell("mkdir -p " + os.path.dirname(rgn) )
			open(rgn, "w+").write(line)	
#
# helper functions for getting inputs/params for SDA local assemblies
#
def get_ids(wildcards):
	checkpoint_output = checkpoints.local_asm_dirs.get(**wildcards).output.LAs
	PRE = wildcards.PRE
	DIR = wildcards.DIR
	TMPS = glob_wildcards( os.path.join(checkpoint_output, "region_{LA_ID}/rgn.bed" ) ).LA_ID
	LA_IDs = []
	# filter for only IDs that are \d+ and convert to ints
	for LA_ID in TMPS:
		if( re.match("\d+", LA_ID)):
			LA_IDs.append(int(LA_ID))
	# sort cuts to garuntee the same ordering each times
	LA_IDs = sorted(LA_IDs)
	# assert that all the cut IDs that should be there are, (no numbers are skipped)
	for idx, val in enumerate(LA_IDs):
		assert idx == val
	return(LA_IDs)

def get_rgn(wildcards):	
	LA_ID = int(str(wildcards.LA_ID))
	token = open(f"{DIR}/coverage/{PRE}.collapses.bed").readlines()[LA_ID].strip().split()
	return("{}:{}-{}".format(token[0], token[1], token[2]))

def get_dir(wildcards):	
	LA_ID = int(str(wildcards.LA_ID))
	return(COL_DIR.format(DIR=DIR, PRE=PRE, LA_ID=LA_ID))

def get_pre(wildcards):	
	LA_ID = int(str(wildcards.LA_ID))
	return(f"region_{LA_ID}")

def get_cov(wildcards):	
	stats = pd.read_csv(f"{DIR}/coverage/{PRE}.coverage.stats", sep="\t")
	mean = stats["mean_coverage"][0]
	return(mean)

#
# rules to gather data for and run SDA
#

rule la_ref:
	input:
		ref = REF,
		hcr = rules.high_coverage_regions.output.hcr,
	output:
		ref = COL_REF_FMT,
		fai = COL_REF_FMT + ".fai",
	params:
		rgn = get_rgn,
	resources:
		mem=16,
	threads:1
	shell:"""
samtools faidx {input.ref} '{params.rgn}' > {output.ref}
samtools faidx {output.ref}
"""

def get_refs(wildcards):
	return( [ COL_REF_FMT.format(DIR=DIR, PRE=PRE, LA_ID=LA_ID) for LA_ID in get_ids(wildcards) ] )

rule la_bam:
	input:
		hcr = rules.high_coverage_regions.output.hcr,
		bam=rules.index_bam.input.bam,
		bai=rules.index_bam.output.bai,
	output:
		bam = COL_BAM_FMT,
	params:
		rgn = get_rgn,
	resources:
		mem=16,
	threads:1
	shell:"""
samtools view -F 2308 -b {input.bam} '{params.rgn}' > {output.bam}
"""

def get_bams(wildcards):
	return( [ COL_BAM_FMT.format(DIR=DIR, PRE=PRE, LA_ID=LA_ID) for LA_ID in get_ids(wildcards) ] )


rule la_sda:
	input:
		bam = COL_BAM_FMT,
		ref = COL_REF_FMT,
		fai = COL_REF_FMT + ".fai",
	output:
		sda = COL_SDA_FMT,
		out = COL_SDA_FMT + ".log",
	params:
		sda_dir = get_dir,
		pre = get_pre,
		cov = get_cov,
	resources:
		mem=4,
	threads:8
	run:
		SDA_BAM = os.path.abspath(input["bam"])
		SDA_REF = os.path.abspath(input["ref"])
		SDA_DIR = os.path.abspath(params["sda_dir"])
		SDA_SDA = os.path.abspath(output["sda"])
		SDA_OUT = os.path.abspath(output["out"])
		cmd = """
# clear any previous runs 
rm -rf {params.sda_dir}/{params.pre}*
rm -rf {params.sda_dir}/*/{params.pre}*

# move to execution dir 
pushd {params.sda_dir}

{snake_dir}SDA collapse --ref {SDA_REF} --reads {SDA_BAM} --coverage {params.cov} \
	-d {SDA_DIR} -p {params.pre} -t {threads} \
	--platform {PLAT} --minaln {MINALN} --bandwidth {BANDWIDTH} --iterations {ITERATIONS} \
	--assemblers {ASSEMBLERS} --lrt {LRT} --minNumShared {MINNUMSHARED} --maxPosRep {MAXPOSREP} \
	--minCutSize {MINCUTSIZE} --minCutLen {MINCUTLEN} --unlock &> \
	/dev/null || echo "Already unlocked."


timeout {MAX_TIME} {snake_dir}SDA collapse --ref {SDA_REF} --reads {SDA_BAM} --coverage {params.cov} \
	-d {SDA_DIR} -p {params.pre} -t {threads} \
	--platform {PLAT} --minaln {MINALN} --bandwidth {BANDWIDTH} --iterations {ITERATIONS} \
	--assemblers {ASSEMBLERS} --lrt {LRT} --minNumShared {MINNUMSHARED} --maxPosRep {MAXPOSREP} \
	--minCutSize {MINCUTSIZE} --minCutLen {MINCUTLEN}  &> \
	{SDA_OUT} || echo "SDA failed on this collapse" && touch {SDA_SDA}

popd
		"""
		shell(cmd)
		

def get_sda(wildcards):
	return( [ COL_SDA_FMT.format(DIR=DIR, PRE=PRE, LA_ID=LA_ID) for LA_ID in get_ids(wildcards) ] )

rule merge_results:
	input:
		get_sda,
	output:
		fasta = "{DIR}/{PRE}.assemblies.fasta",
		preads = "{DIR}/{PRE}.phased.readids",
		summary = "{DIR}/{PRE}.summary.txt",
		psvs = "{DIR}/{PRE}.psv.tbl",
	resources:
		mem=4,
	threads: 1
	run:
		preads = []
		summary = []
		for LA_ID, sda_log in enumerate(input):	
			print(LA_ID, sda_log)
			sda_dir = os.path.dirname(sda_log)
			
			fasta_path = os.path.join(sda_dir, f"region_{LA_ID}.assemblies.fasta")
			if(os.path.exists(fasta_path)):
				shell("cat {fasta_path} >> {output.fasta}")

			preads_path = os.path.join(sda_dir, f"region_{LA_ID}.phased.readids")
			if(os.path.exists(preads_path)):
				preads.append(pd.read_csv(preads_path, sep = "\t"))

			summary_path = os.path.join(sda_dir, f"region_{LA_ID}.summary.txt")
			if(os.path.exists(summary_path)):
				summary.append(pd.read_csv(summary_path, sep = "\t"))

			psv_path = os.path.join(sda_dir, f"region_{LA_ID}.psv.tbl")
			if(os.path.exists(psv_path)):
				shell("cat {psv_path} >> {output.psvs}")
		
		pd.concat(preads, ignore_index=True).to_csv(output["preads"], sep="\t", index=False)
		pd.concat(summary, ignore_index=True).to_csv(output["summary"], sep="\t", index=False)

rule summary_plots:
	input:
		fasta = "{DIR}/{PRE}.assemblies.fasta",
		preads = "{DIR}/{PRE}.phased.readids",
		summary = "{DIR}/{PRE}.summary.txt",
		psvs = "{DIR}/{PRE}.psv.tbl",
	output:
		pair = "{DIR}/summary_plots/{PRE}.paired.pdf",
		length = "{DIR}/summary_plots/{PRE}.assembled_lengths.pdf",
		collapse = "{DIR}/summary_plots/{PRE}.collapse_vs_assemblies.pdf",
		bar = "{DIR}/summary_plots/{PRE}.assembled_mbp.pdf",
	resources:
		mem=16, 
	threads: 1
	shell:"""
{snake_dir}scripts/SummarySDAPlots.py  {input.summary} \
	--dir {DIR} \
	--prefix {PRE} \
	--pair {output.pair} \
	--length {output.length} \
	--collapse {output.collapse} \
	--bar {output.bar}
"""

rule final:
	input:
		fasta = "{DIR}/{PRE}.assemblies.fasta",
		preads = "{DIR}/{PRE}.phased.readids",
		summary = "{DIR}/{PRE}.summary.txt",
		psvs = "{DIR}/{PRE}.psv.tbl",
		pair = rules.summary_plots.output.pair,
		length = rules.summary_plots.output.length,
		collapse = rules.summary_plots.output.collapse,
		bar = rules.summary_plots.output.bar,
	output:
		final = "{DIR}/{PRE}.done",
	resources:
		mem=16,
	threads:1
	shell:"""
touch {output.final}
"""

 
