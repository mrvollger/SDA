import os
import sys
import tempfile 
import numpy as np
import re
import pysam
from collections import Counter
from Bio import SeqIO
import pandas as pd

snake_dir = os.path.dirname(workflow.snakefile) + "/"
CWD = os.getcwd()

python2 = snake_dir + "envs/env_python2.sh" 
python3 = snake_dir + "envs/env_python3.sh" 

shell.executable("/bin/bash")
shell.prefix(f"cd {snake_dir}/envs && source {python3} && cd {CWD} && set -eo pipefail; ")

base = snake_dir + "scripts/"
scriptsDir = snake_dir + "CCscripts/"
CANU = snake_dir + "externalRepos/canu-1.8/Linux-amd64/bin/canu"
RACON =snake_dir + "externalRepos/racon-v1.4.5/build/bin/racon"

# script locations and configurations 
#
COV = config["coverage"]
PLAT = config["platform"].upper()
READS = os.path.abspath( config["reads"] )
REF = os.path.abspath( config["ref"] )
MINALN = config["minaln"] 
BWIDTH = config["bandwidth"] 
PRE	= config["prefix"]
DIR	= config["dir"]
MINLRT = config["lrt"]
MINSHARED = config["minNumShared"]
MAXPOSREP = config["maxPosRep"]
MINCUTSIZE = config["minCutSize"]
MINCUTLEN = config["minCutLen"]
ITERATIONS = config["iterations"]
DEBUG=config["debug"]
MINASMLENGTH=10000

# standard deviation of read depth
STD = np.sqrt(COV)
MINCOV = int( max(COV - 4*STD, COV/2.0) )
MAXCOV = int( COV + 1*STD )
MINTOTAL = int( 2*COV - 3*STD )
MINREADS = int(MINCOV/2.0)


ASMS = config["assemblers"].strip().split(",")

wildcard_constraints:
	PRE = PRE,
	DIR = DIR,
	CUT = "\d+",
	ASM = "|".join(ASMS)

assert os.path.exists(REF+".fai"), "{} must be indexed!".format(REF)

onsuccess:
	sys.stderr.write("SDA FINISHED\n")
onerror:
	sys.stderr.write("SDA FAILED\n")

rule all:	
	input:
		done=expand("{DIR}/{PRE}.done",DIR=DIR, PRE=PRE),
		
def tempd(File):
	if(DEBUG):
		return(File)
	return(temp(File))



if(PLAT in ["CCS", "SUBREAD"] ):
	rule pbmm2:
		input:
			reads = READS,
			ref = REF,
		output:
			tempd("{DIR}/{PRE}.tmp.reads.bam")
		threads: 8   
		shell: """ 
# SUBREAD gives bettern alignmetns for CCS reads than CCS
pbmm2 align -j {threads} \
	--preset SUBREAD -N 50  --min-length {MINALN} -r {BWIDTH} \
	--sample FAKE_SAMPLE \
	{input.ref} {input.reads} | \
	samtools view -u -F 2308 - | \
	samtools sort -@ {threads} -m 4G -T tmp -o {output}
"""

elif(PLAT in ["ONT"] ): 
	rule minimap2:
		input:
			reads = READS,
			ref = REF,
		output:
			tempd("{DIR}/{PRE}.tmp.reads.bam")
		threads: 8
		shell:"""
samtools fasta {input.reads} | \
	minimap2 \
	-ax map-ont \
	--eqx -L \
	-R '@RG\\tID:MINIMAP\\tSM:FAKE_SAMPLE\\tPL:ONT' \
	-t {threads} \
	-m {MINALN} -r {BWIDTH} \
	{input.ref} /dev/stdin | \
	samtools view -u -F 2308 - | \
	samtools sort -@ {threads} -m 4G -T tmp -o {output}
"""
else:
	sys.stderr.write("Platform {} not recongnized!\n".format(PLAT))


# Marks readToSNVList requires that tlen in the bam be set to the read reference length 
# so I must convert the bam 
rule change_bam:
	input:
		bam="{DIR}/{PRE}.tmp.reads.bam"
	output:
		bam = tempd("{DIR}/{PRE}.reads.bam"),
	shell:"""
{base}changeBamName.py {input.bam} {output.bam}
"""

rule index:
	input:
		bam = "{DIR}/{PRE}.reads.bam",
	output:
		bai = tempd("{DIR}/{PRE}.reads.bam.bai"),
	shell:"""
samtools index {input.bam}
"""

rule nucplot:
	input:
		bam = rules.index.input.bam,
		bai = rules.index.output.bai,
	output:
		png = "{DIR}/{PRE}.coverage.png",
	shell:"""
{snake_dir}externalRepos/nucfreq/NucPlot.py {input.bam} {output.png}
"""
		

#
# Given the alignments, count the frequency of each base at every
# position, and estimate the SNVs from this frequency. 
#
rule nucfreq:
	input:
		bam = rules.index.input.bam,
		bai = rules.index.output.bai,
		ref= REF,
	output:
		nucfreq=tempd("{DIR}/snvs/{PRE}.assembly.consensus.nucfreq"),
	shell:"""
source {python2}
echo "Sam to nucfreq"
samtools mpileup -q 0 -Q 0 {input.bam} | \
    {base}/MpileupToFreq.py  /dev/stdin | \
    {base}/PrintHetFreq.py {MINCOV} \
    --maxCount {MINCOV} \
    --minTotal {MINTOTAL} \
    > {output.nucfreq}

if [ -s {output.nucfreq} ]; then
    echo "Candidate PSVs found!"
else
    >&2 echo "ERROR: NO candidate PSVS were found. SDA cannot resolve this duplication."
    rm {output}
    exit 1
fi
"""


rule snv_tbl:
	input:
		bam = rules.index.input.bam,
		nucfreq=rules.nucfreq.output.nucfreq,
		ref = REF,
	output: 
		snv= tempd("{DIR}/snvs/{PRE}.assembly.consensus.fragments.snv"),
		vcf= tempd("{DIR}/snvs/{PRE}.assembly.consensus.nucfreq.vcf"),
		filt=tempd("{DIR}/snvs/{PRE}.assembly.consensus.nucfreq.filt"),
		frag=tempd("{DIR}/snvs/{PRE}.assembly.consensus.fragments"),
	shell:"""
source {python2}
echo "filter nucfreq"
# this actaully only really creates the frag file, all filtering should be done by PrintHetfreq
samtools view -h {input.bam} | \
	{base}/readToSNVList  \
	--nft {input.nucfreq} \
	--sam /dev/stdin \
	--ref {input.ref} \
	--minFraction 0.01 \
	--minCoverage {MINCOV} \
	--out {output.frag} \
	--nftOut {output.filt}

echo "filtered to vcf"
{scriptsDir}/FreqToSimpleVCF.py --freq {output.filt} \
	--ref {input.ref} \
	--out {output.vcf} 

echo "vcf to snv"
{scriptsDir}/FragmentsToSNVList.py  \
	--fragment \
	--frags {output.frag} \
	--vcf {output.vcf} \
	--out {output.snv}
"""



#
# Create a matrix with one row per read, and one column per PSV.  
#  . - read is ref base
#  1 - read is PSV
#  n - no signal (indel or read ended)
#
rule snv_mat:
    input:
        snv=rules.snv_tbl.output.snv,	
    output:
        mattmp=tempd("{DIR}/snvs/{PRE}.assembly.consensus.fragments.snv.mat"),
        mat= tempd("{DIR}/snvs/{PRE}.assembly.consensus.fragments.snv.mat.catergorized"),
        snvpos=tempd("{DIR}/snvs/{PRE}.assembly.consensus.fragments.snv.pos")
    shell:"""
source {python2}
{scriptsDir}/FragmentSNVListToMatrix.py {input.snv} --named --pos {output.snvpos} --mat {output.mattmp}

# this adds a fake catagory on the end, becasue downstream scripts need that column to exist
cat {output.mattmp} | awk '{{ print $1"\t"$2"\tall"}}' > {output.mat}
"""


#
# This finds PSVs that are connected by a sufficient number of
# sequences, and creates the PSV graph. This will have merged components.
#
rule psv_graph:
	input:
		mat=rules.snv_mat.output.mat,
		vcf=rules.snv_tbl.output.vcf, 
	output:
		graph=tempd("{DIR}/CC/{PRE}.mi.gml"),
		adj = tempd("{DIR}/CC/{PRE}.mi.adj"),
		mi  = tempd("{DIR}/CC/{PRE}.mi.mi"),
	shell:"""
source {python2}
{scriptsDir}/PairedSNVs.py {input.mat} --minCov {MINCOV} --maxCov {MAXCOV} \
		--mi {output.mi} --graph {output.graph} --adj {output.adj} \
		--minNShared {MINSHARED} --minLRT {MINLRT} --vcf {input.vcf}

"""



#
# generate repulsion edges 
#
rule repulsion:
	input:
		graph=rules.psv_graph.output.graph,
		mi= rules.psv_graph.output.mi,
	output:
		rep = tempd("{DIR}/CC/{PRE}.mi.repulsion"),
	shell:"""
source {python2}
{base}/GenerateRepulsion.py --shared {MINSHARED} --lrt {MINLRT} --max {MAXPOSREP} --gml {input.graph} --mi {input.mi} --out {output.rep}
"""

#
# Run correlation clustering to try and spearate out the graph.  Swap
# is a 'simulated annealing' type parameters. The factor parameter
# controls for how many negative edges are added for every positive edge.
#
rule run_cc:
	input:
		graph=rules.psv_graph.output.graph,
		rep=rules.repulsion.output.rep, 
	output:
		graph="{DIR}/CC/{PRE}.mi.cuts.gml",
		sites="{DIR}/CC/{PRE}.mi.gml.sites",
		cuts="{DIR}/CC/{PRE}.mi.gml.cuts",
		score="{DIR}/CC/{PRE}.cuts.score.txt",
	threads: 8
	shell:"""	
source {python2}

{scriptsDir}/MinDisagreeClusterByComponent.py  \
	--plotRepulsion \
	--graph {input.graph} \
	--repulsion {input.rep} \
	--niter 10000 --swap 1000000 --factor 1 \
	--minCutSize {MINCUTSIZE} \
	--minlen {MINCUTLEN} \
	--starts {ITERATIONS} \
	--threads {threads} \
	--scores {output.score} \
	--cuts {output.cuts} --sites {output.sites} --out {output.graph}

if [ -s {output.sites} ]; then
	echo "CC groups created!"
else
	>&2 echo "ERROR: No CC groups were made. SDA found no clusters of PSVs and cannot resolve this duplication."
	rm {output}
	exit 1
fi
"""

#
# makes a modified coverage plot with CC groups drawn on top
#
rule nucplot_cc:
	input:
		bam = rules.index.input.bam,
		bai = rules.index.output.bai,
		sites=rules.run_cc.output.sites, 
	output: 
		png="{DIR}/{PRE}.CovWithCC.png",
	shell:"""
{snake_dir}externalRepos/nucfreq/NucPlot.py --psvsites {input.sites} {input.bam} {output.png}
"""


#
# makes a gephi version of the plot, 
#
rule gephi:
	input:
		graph=rules.run_cc.output.graph, 
	output:
		pdf="{DIR}/CC_plots/{PRE}.cc.pdf",
	shell:"""
{base}/gephi/gephi.sh {input.graph} {output.pdf}
"""


#
# Correlation clustering defines a set of cuts that separate out
# connected components.  This takes the cuts and makes a vcf file for
# each component.
#
checkpoint cuts_to_vcfs:
	input:
		refIdx=REF + ".fai",
		cuts=rules.run_cc.output.cuts, 
		snvpos=rules.snv_mat.output.snvpos,
		vcf=rules.snv_tbl.output.vcf, 
	output:
		#vcfs=dynamic("group.{n}.vcf"), # this caused a snakemake effore often because dynamic is depricated 
		#done = "CC/vcfs.done",
		groups = directory("{DIR}/{PRE}.cuts/"),
		comps="{DIR}/CC/{PRE}.mi.comps.txt", 
	shell:"""
source {python2}
{scriptsDir}/CutsToPhasedVCF.py {input.cuts} {input.snvpos} {input.vcf} \
		--minComponent {MINCUTSIZE} \
		--ref {input.refIdx} \
		--base {output.groups}/{wildcards.PRE} \
		--summary {output.comps}
"""

def get_cuts(wildcards):
	checkpoint_output = checkpoints.cuts_to_vcfs.get(**wildcards).output.groups
	PRE = wildcards.PRE
	DIR = wildcards.DIR
	TMPS = glob_wildcards( os.path.join(checkpoint_output, "{PRE}.{{CUT}}.vcf".format( PRE=PRE ) )   ).CUT
	CUTS = []
	# filter for only cuts that are \d+ and convert to ints
	for cut in TMPS:
		if( re.match("\d+", cut)):
			CUTS.append(int(cut))
	# sort cuts to garuntee the same ordering each times
	CUTS = sorted(CUTS)
	# assert that all the cut IDs that should be there are, (no numbers are skipped)
	for idx, val in enumerate(CUTS):
		assert idx == val 
	return(CUTS)



#-----------------------------------------------------------------------------------------------------#
#
# this part partitions the reads based on the vcf files (PSVs)
#
# make a phased vcf for whatshap, just a format change
rule phased_vcf:
	input: 
		vcf = '{DIR}/{PRE}.cuts/{PRE}.{CUT}.vcf',
		bam = rules.index.input.bam,
		bai = rules.index.output.bai,
	output: 
		vcf = tempd("{DIR}/{PRE}.cuts/{PRE}.phased.{CUT}.vcf")
	shell:"""
{base}/fixVCF.py --out {output.vcf} --vcf {input.vcf} --bam {input.bam}
"""

# run whats hap and get the partitioned sam files
rule whatshap:
	input: 
		bam = rules.index.input.bam,
		bai = rules.index.output.bai,
		vcf = rules.phased_vcf.output.vcf,
	output: 
		bam= '{DIR}/{PRE}.cuts/{PRE}.{CUT}.H2.WH.bam'
	shell:"""
whatshap haplotag {input.vcf} {input.bam} | \
	samtools view -h - | \
	grep -E "^@|HP:i:2" | \
	samtools view -b - > {output.bam}
"""



#
#-----------------------------------------------------------------------------------------------------#
#
# This part runs the assembly based on the partitions 
#

# get fasta files from the reads
rule phased_reads:
	input: 
		bam = rules.whatshap.output.bam, 
	output:
	    reads=tempd('{DIR}/{PRE}.cuts/{PRE}.{CUT}.reads.fasta')
	shell:"""
samtools fasta {input.bam} > {output.reads} 
"""


# check the input sequence data for the purpose of canu 
if(PLAT == "CCS"):
	canu_p = f' corMinCoverage=0 contigFilter="{MINREADS} {MINASMLENGTH} 1.0 .25 {MINREADS}" ovlMerThreshold=75 batOptions="-eg 0.01 -eM 0.01 -dg 6 -db 6 -dr 1 -ca 50 -cp 5"  correctedErrorRate=0.015 -pacbio-corrected '
	wtdbg2_p = " ccs "
elif(PLAT == "SUBREAD"):
	canu_p = f' corMinCoverage=1 contigFilter="{MINREADS} {MINASMLENGTH} 1.0 .75 {MINREADS}" -pacbio-raw '
	wtdbg2_p = " sequel "
elif(PLAT=="ONT"):
	canu_p = f' corMinCoverage=1 contigFilter="{MINREADS} {MINASMLENGTH} 1.0 .75 {MINREADS}" -nanopore-raw '
	wtdbg2_p = " ont "


# run the assembly
rule run_asm:
	input: 
		rules.phased_reads.output.reads, 
	output: 
		asm = tempd('{DIR}/{PRE}.cuts/{PRE}.{CUT}.{ASM}.fasta'),
	threads: 4
	run:
		ASM = str(wildcards.ASM)
		DIR = str(wildcards.DIR)
		PRE = str(wildcards.PRE)
		CUT = int(str(wildcards.CUT))
		# asm location
		asmdir=f'{DIR}/{PRE}.cuts/{PRE}.{CUT}.{ASM}/'
		print(asmdir)
		shell("rm -rf {asmdir} && mkdir -p {asmdir}")
	
		if(ASM == "canu"):
			shell('''
				{CANU} {canu_p} {input} \
						stopOnLowCoverage=1 \
						genomeSize=60000 \
						corOutCoverage=300 \
						corMhapSensitivity=high \
						-p asm useGrid=false  \
						-d {asmdir} \
						maxThreads={threads} \
						&& mv {asmdir}/asm.contigs.fasta {output} \
						|| ( >&2 echo " no real assembly" && > {output} )
				''')
		elif(ASM == "wtdbg2"):
			shell('''
				wtdbg2 -f -i {input} \
					-x {wtdbg2_p} \
					-o {asmdir}/asm  \
					--ctg-min-length {MINASMLENGTH} \
					-t {threads} \
					-L 5000  

				# run consensous
				if [ -s {asmdir}/asm.ctg.lay.gz ]; then
					wtpoa-cns \
							-t {threads} \
							-i {asmdir}/asm.ctg.lay.gz \
							-fo {output}
				else
					touch {output}
				fi
				''')
			
		shell("rm -rf {asmdir} ")


def get_asms(wildcards):
	return( expand("{DIR}/{PRE}.cuts/{PRE}.{CUT}.{ASM}.fasta", DIR=wildcards.DIR, CUT=get_cuts(wildcards), PRE=wildcards.PRE, ASM=ASMS ) )

#-----------------------------------------------------------------------------------------------------#
#
# run polishing steps
#

rule prep_for_polish:
	input:
		asm = rules.run_asm.output.asm, 
		bam = rules.whatshap.output.bam, 
	output: 
		bam = tempd('{DIR}/{PRE}.cuts/{PRE}.{CUT}.{ASM}.bam'),
		bai = tempd('{DIR}/{PRE}.cuts/{PRE}.{CUT}.{ASM}.bam.bai'),
		pbi = tempd('{DIR}/{PRE}.cuts/{PRE}.{CUT}.{ASM}.bam.pbi'),
	threads: 4
	run:
		if(os.path.getsize( input["asm"] ) == 0):
			shell("touch {output}")	
		elif(PLAT in ["CCS", "SUBREAD"]):
			shell("pbmm2 align --preset SUBREAD -j {threads} {input.asm} {input.bam} | samtools view -u -F 2308 - | samtools sort - > {output.bam}" )
			shell("samtools index {output.bam} && pbindex {output.bam}")
		else:
			shell("samtools fastq {input.bam} | minimap2 -ax map-ont --eqx -L -t {threads} -m {MINALN} -r {BWIDTH} {input.asm} /dev/stdin | samtools view -u -F 2308 - | samtools sort - > {output.bam}" )
			shell("samtools index {output.bam} && touch {output.pbi}")


rule polish_asm:
	input:
		asm = rules.run_asm.output.asm, 
		bam = rules.prep_for_polish.output.bam,
		bai = rules.prep_for_polish.output.bai,
		pbi = rules.prep_for_polish.output.pbi,
	output:
		fasta = tempd('{DIR}/{PRE}.cuts/{PRE}.{CUT}.{ASM}.cor.fasta'),
		tmp_sam = tempd('{DIR}/{PRE}.cuts/{PRE}.{CUT}.{ASM}.sam'),
		tmp_fastq = tempd('{DIR}/{PRE}.cuts/{PRE}.{CUT}.{ASM}.reads.fastq'),
		tmp_cor = tempd('{DIR}/{PRE}.cuts/{PRE}.{CUT}.{ASM}.racon.fasta'),
		tmp_fai = tempd('{DIR}/{PRE}.cuts/{PRE}.{CUT}.{ASM}.racon.fasta.fai'),
	threads: 4
	run:
		if(os.path.getsize( input["asm"] ) == 0):
			shell("touch {output}")	
		else:
			# make inputs for racon
			shell("samtools view {input.bam} > {output.tmp_sam}")
			shell("samtools fastq {input.bam} > {output.tmp_fastq}")
			
			# run racon ( or try to )
			shell("""{RACON} {output.tmp_fastq} {output.tmp_sam} {input.asm} -u -t {threads} > {output.tmp_cor} || \
				( >&2 echo " RACON FAILED TO RUN " && cp {input.asm} {output.tmp_cor} ) """)	

			# run arrow ( or try to ) 
			shell("samtools faidx {output.tmp_cor}")	
			if(PLAT in ["SUBREAD"]):
				shell("""source {python2} && \
					arrow --noEvidenceConsensusCall lowercasereference \
					--minCoverage 10 -j {threads} -r {output.tmp_cor} -o {output.fasta} {input.bam} || \
					( >&2 echo " ARROW FAILED TO RUN " && cp {output.tmp_cor} {output.fasta} )
					""")
			else:
				shell("cp {output.tmp_cor} {output.fasta}")

def get_polish(wildcards):
	return( expand("{DIR}/{PRE}.cuts/{PRE}.{CUT}.{ASM}.cor.fasta", DIR=wildcards.DIR, CUT=get_cuts(wildcards), PRE=wildcards.PRE, ASM=ASMS ) )


#-----------------------------------------------------------------------------------------------------#


def write_fastas(d_recs, outfile, counter = 0):
		out = []
		for cut, asm in d_recs:
			recs = d_recs[(cut,asm)]
			for rec in recs:
				rec.id = f"prefix.{PRE}_group.{cut}_asm.{asm}_uid.{counter}"
				rec.name = rec.id
				counter += 1	
				out.append(rec)
		SeqIO.write(out, outfile, "fasta")
		return(counter)

rule merge_asms:
	input:
		asms = get_polish,
	output:
		asms = "{DIR}/{PRE}.assemblies.fasta",
		unused = "{DIR}/{PRE}.unused.assemblies.fasta",
	run:
		PRE = str(wildcards["PRE"])
		asms = "|".join(ASMS)
		recs = {}
		for asm in ASMS:
			recs[asm] = {}
		
		pat = f'{DIR}/{PRE}\.cuts/{PRE}\.(\d+)\.({asms})\.cor\.fasta'
		for fasta in input["asms"]:
			match = re.match(pat, fasta)
			assert match 
			cut, asm = match.groups()
			cut = int(cut)
			recs[asm][cut] = list( SeqIO.parse(fasta, "fasta") )	
		
		rtn = {}
		unused = {}
		for cut in get_cuts(wildcards):
			saved = False
			for asm in ASMS:	
				fasta = recs[asm][cut]
				if( len(fasta) > 0 and not saved):
					rtn[(cut, asm)] = fasta
					saved = True
				else:
					unused[(cut, asm)] = fasta

		counter = write_fastas(rtn, output["asms"])
		counter = write_fastas(unused, output["unused"], counter=counter)
		

rule aln_asms:
	input:
		ref = REF,
		asms = rules.merge_asms.output.asms,
	output:
		bam = tempd("{DIR}/{PRE}.asm.bam"),
		bai = tempd("{DIR}/{PRE}.asm.bam.bai"),
	threads: 4
	shell:"""
minimap2 -t {threads} --secondary=no -a --eqx -Y -x asm20 \
	-m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5 \
	{input.ref} {input.asms} | samtools view -F 2308 -u - | samtools sort -@ {threads} - > {output.bam}
samtools index {output.bam}
"""


def get_vcfs(wildcards):
	return( expand("{DIR}/{PRE}.cuts/{PRE}.{CUT}.vcf", DIR=wildcards.DIR, CUT=get_cuts(wildcards), PRE=wildcards.PRE) )

rule psv_pos:
	input:
		vcfs=get_vcfs,
		asms= rules.merge_asms.output.asms,
		bam = rules.aln_asms.output.bam,
		bai = rules.aln_asms.output.bai,
	output:
		psvtbl = "{DIR}/{PRE}.psv.tbl",
	shell:"""
{base}PSVLocations.py --check {input.asms} --bam {input.bam} --psvs {input.vcfs} > {output.psvtbl}
"""




def get_phased_bams(wildcards):
	return( expand("{DIR}/{PRE}.cuts/{PRE}.{CUT}.H2.WH.bam", DIR=wildcards.DIR, CUT=get_cuts(wildcards), PRE=wildcards.PRE ) )

rule reads_by_cut:
	input:
		bams = get_phased_bams,
	output:
		pids = "{DIR}/{PRE}.phased.readids",
	run:
		names = []
		nameByCut = {}
		
		# read in basm enumerated by cuts 	
		for cut, bam in enumerate(input["bams"]):
			bam = pysam.AlignmentFile(bam)
			for idx, read in enumerate(bam.fetch(until_eof=True)):
				if(cut not in nameByCut):
					nameByCut[cut] = set()
				nameByCut[cut].add(read.query_name)
				names.append(read.query_name)
			bam.close()

		counts = Counter(names)
		out = "prefix\treadid\tcut\tunique\n"
		for cut in sorted(nameByCut):
			for readid in nameByCut[cut]:
				unique = True
				if(counts[readid] > 1):
					unique = False
				out += f"{PRE}\t{readid}\t{cut}\t{unique}\n"
		
		open(output["pids"], "w+").write(out)


rule summary:
	input:
		pids = rules.reads_by_cut.output.pids,
		sites=rules.run_cc.output.sites, 
		asms = rules.merge_asms.output.asms,
	output:
		summary = "{DIR}/{PRE}.summary.txt",
	run:
		pids = pd.read_csv(input["pids"], sep="\t")	
		psvs = {"cut":[], "nPSVs":[], "nReads":[], "nUniqueReads":[]}
		for cut, line in enumerate(open(input["sites"])):
			PSVs = line.strip().split()
			nPSVs=len(PSVs)
			nReads = sum( pids["cut"] == cut)
			nUniqueReads = sum( (pids["cut"] == cut) & ( pids["unique"] ) )	
			psvs["cut"].append(cut)
			psvs["nPSVs"].append(nPSVs)
			psvs["nReads"].append(nReads)
			psvs["nUniqueReads"].append(nUniqueReads)
		psvs = pd.DataFrame(psvs)	
		
		pasms = "|".join(ASMS)	
		pat = f'prefix\.{PRE}_group\.(\d+)_asm\.({pasms})_uid\.(\d+)'
		asms = {"uid":[], "cut":[], "length":[], "assembler":[] , "asm_name":[]}
		for rec in SeqIO.parse(input["asms"], "fasta"):
			match = re.match(pat, rec.id)
			assert match, (pat, rec.id)
			cut, asm, uid = match.groups()
			asms["uid"].append(int(uid))
			asms["cut"].append(int(cut))
			asms["length"].append( len(rec.seq) )	
			asms["assembler"].append( asm )	
			asms["asm_name"].append( rec.id )	
		asms = pd.DataFrame(asms)
	
		rtn = pd.merge(psvs, asms, on='cut', how='outer') 
		rtn["asm_failed"] = pd.isna( rtn["length"] )
		rtn["multiple_asm"] =  rtn["cut"].duplicated(keep=False)
		
		# reorder the columns 
		cols = list(rtn)
		cols.insert(1, cols.pop(cols.index("uid")))
		cols.insert(len(cols), cols.pop(cols.index("asm_name")))
		rtn = rtn[cols]
		
		# add prefix
		rtn["prefix"] =  f"{PRE}"

		pd.set_option('precision', 0); print(rtn)
		rtn.to_csv(output["summary"], sep="\t", index=False)




#-----------------------------------------------------------------------------------------------------#
#
# Plot phased reads agaisnt asm
#
rule prep_for_plot:
	input:
		asm = rules.polish_asm.output.fasta, 
		bam = rules.whatshap.output.bam, 
	output: 
		bam = tempd('{DIR}/{PRE}.cuts/{PRE}.{CUT}.cor_{ASM}.bam'),
		bai = tempd('{DIR}/{PRE}.cuts/{PRE}.{CUT}.cor_{ASM}.bam.bai'),
	threads: 4
	run:
		if(os.path.getsize( input["asm"] ) == 0):
			shell("touch {output}")	
		elif(PLAT in ["CCS", "SUBREAD"]):
			shell("pbmm2 align --preset SUBREAD -j {threads} {input.asm} {input.bam} | samtools view -u -F 2308 - | samtools sort - > {output.bam}" )
			shell("samtools index {output.bam}")
		else:
			shell("samtools fastq {input.bam} | minimap2 -ax map-ont --eqx -L -t {threads} -m {MINALN} -r {BWIDTH} {input.asm} /dev/stdin | samtools view -u -F 2308 - | samtools sort - > {output.bam}" )
			shell("samtools index {output.bam}")

rule plot_asm:
	input:
		bam = rules.prep_for_plot.output.bam,
		bai = rules.prep_for_plot.output.bai,
	output: 
		png = "{DIR}/asm_plots/{PRE}.{CUT}.{ASM}.png"	
	threads: 4
	run:
		if(os.path.getsize( input["bam"] ) == 0):
			shell("touch {output}")	
		else:
			shell("{snake_dir}externalRepos/nucfreq/NucPlot.py {input.bam} {output.png} --ylim {MINTOTAL}")


def get_phased_plots(wildcards):
	return( expand("{DIR}/asm_plots/{PRE}.{CUT}.{ASM}.png", DIR=wildcards.DIR, CUT=get_cuts(wildcards), PRE=wildcards.PRE, ASM=ASMS ) )


rule final_rule:
	input:
		plots = get_phased_plots,
		asms = rules.merge_asms.output.asms,
		psvs = rules.psv_pos.output.psvtbl, 
		summary = rules.summary.output.summary,
		png = rules.nucplot.output.png, 	
		cuts = rules.run_cc.output.cuts,	
		png2 = rules.nucplot_cc.output.png, 
		pdf = rules.gephi.output.pdf,
		comps = "{DIR}/CC/{PRE}.mi.comps.txt",
	output:
		done="{DIR}/{PRE}.done",
	shell:"""
touch {output.done}
"""



