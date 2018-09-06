import os
import glob
import re
from Bio import SeqIO
import re


snake_dir = os.path.dirname(workflow.snakefile)+"/"
shell.executable("/bin/bash")
shell.prefix("source {}/env_python3.cfg; ".format(snake_dir))

#
# script locations and configurations 
#
base = snake_dir + "scripts/"
python3 = snake_dir + "env_python3.cfg"
python2 = snake_dir + "env_python2.cfg"


if(os.path.exists("abp.config.json")):
	cfile = "abp.config.json"
else:
	cfile = "sda.config.json"
configfile:
	cfile



ISONT=False
MM2_TYPE=" map-pb " 
if("ont" in config):
	if(config["ont"].lower() in ["t", "true"] ):
		ISONT=True
		MM2_TYPE= " map-ont "
MM2 = False
if("minimap" in config):
	if(config["minimap"].lower() in ["t", "true"] ):
		MM2=True
bandwidth = " 70000 "

# min cov is used to detemrine the filter for getting rid of low read assemblies. 
MINREADS = int(config["MINCOV"]*1.0/2.0)

# learn how many CC clusters there are
groups= glob.glob("group.[0-9]*.vcf")
IDS= []
for group in groups:
    group= group.strip().split(".")
    assert group[0]=="group"
    IDS.append(group[1])


assemblers = ["canu", "wtdbg", "miniasm"]


rule all:
	input: "final",
	message: "Running SDA.2"



# global wild card constraint on n whihc is the group idenitifier
wildcard_constraints:
    n="\d+"


#-----------------------------------------------------------------------------------------------------#
#
# make the group directories, and create a file that jsut acts as a tag for when the dir was created
#
rule makeGroupDirs:
	input:
		vcfs=expand('group.{n}.vcf', n=IDS),
	output: 
		tag=expand('group.{n}/group.{n}', n=IDS)
	run:
        # remove final output
		shell("rm -rf final *.assemblies.fasta")

		for vcf, tag in zip(input["vcfs"], output["tag"] ):
			group, ID, vcf = vcf.split(".")
			curdir = "group.{}/".format(ID)
			shell("rm -rf {}; mkdir -p {}; touch {}".format(curdir, curdir, tag))

#-----------------------------------------------------------------------------------------------------#




#-----------------------------------------------------------------------------------------------------#
#
# this part partitions the reads based on the vcf files (PSVs)
#

# make a phased vcf for whatshap, just a format change
# also requires that pysam is installed, should be taken care of by loading the whatshap anaconda env 
rule phasedVCF:
    input: 'group.{n}/group.{n}',
        vcf= 'group.{n}.vcf',
        hapbam= 'reads.sample.bam',
    output: 'group.{n}/phased.{n}.vcf'
    shell:
        """
		{base}/fixVCF.py --out {output} --vcf {input.vcf} --bam {input.hapbam}
        """

# run whats hap and get the partitioned sam files
rule whatsHap:
    input: 
        hapbam= 'reads.sample.bam',
        hapbai= 'reads.sample.bam.bai',
        hapvcf= 'group.{n}/phased.{n}.vcf'
    output: 
        #hap= 'group.{n}/haplotagged.bam',
        #hapH1= 'group.{n}/H1.WH.sam',
        hapH2= 'group.{n}/H2.WH.bam'
    shell:
        """
		whatshap haplotag {input.hapvcf} {input.hapbam} | \
				samtools view -h - | \
				grep -E "^@|HP:i:2" | \
				samtools view -bS - > {output.hapH2}
        """
#-----------------------------------------------------------------------------------------------------#

rule duplication_index:
	input: 
		expand("group.{n}/H2.WH.bam", n=IDS),
	output:
		dupindex="read.duplication.index",
	shell:"""
{base}/ReadDuplicationIndex.py {input} --out {output.dupindex}
"""


#
#-----------------------------------------------------------------------------------------------------#
#
# This part runs the assembly based on the partitions 
#

# get fasta files from the reads
# if there are no reads partitioned this makes an empty file that subsequent steps check for
rule readsFromSam:
    input: 
        H2= 'group.{n}/H2.WH.bam',
    output:
        pfasta= 'group.{n}/WH.reads.fasta'
    shell:
        """
		grep -v "^@" {input.H2} > group.{wildcards.n}/WH.temp.txt \
			|| touch group.{wildcards.n}/WH.temp.txt

		if [ -s group.{wildcards.n}/WH.temp.txt ]; then
			samtools fasta {input.H2} > {output.pfasta} 
        else
            >&2 echo " no real assembly"
            touch {output.pfasta};
        fi

        rm -f group.{wildcards.n}/WH.temp.txt 
        """

# check the input sequence data for the purpose of canu 
datatype = " -pacbio-raw "
overlaper = " ava-pb "
if(ISONT):
	datatype = " -nanopore-raw "
	overlaper = " ava-ont "

# run the assembly
rule runAssembly:
    input: 'group.{n}/WH.reads.fasta'
    output: 'group.{n}/{ASM}.assembly/asm.contigs.fasta'
    threads: 4
	shell:'''
# asm location
PREFIX="group.{wildcards.n}/{wildcards.ASM}.assembly" 
echo $PREFIX
# make sure we actaully re run the assembly
rm -rf $PREFIX/*
# make asm dir
mkdir -p $PREFIX 
# set minimum assembly length 
MINLENGTH="10000"

if [ -s {input} ]; then
	
	########################## CANU ##########################################
	if [ {wildcards.ASM} == 'canu' ]; then 
		canu {datatype} {input} \
			genomeSize=60000 \
			corOutCoverage=300 \
			corMhapSensitivity=high \
			corMinCoverage=1 \
			gnuplotTested=true  \
			-p asm useGrid=false  \
			-d $PREFIX \
			maxThreads={threads} cnsThreads={threads} ovlThreads={threads} \
			mhapThreads={threads} \
			contigFilter="{MINREADS} $MINLENGTH 1.0 .75 {MINREADS}" \
			|| ( >&2 echo " no real assembly" && \
			> {output} )
	
	########################## WTDBG ##########################################
	elif [ {wildcards.ASM} == 'wtdbg' ]; then 
		wtdbg -i {input} \
				-t {threads} \
				--ctg-min-length $MINLENGTH \
				-o $PREFIX/asm 
		# run consensous
		if [ -s $PREFIX/asm.ctg.lay ]; then
			wtdbg-cns \
					-t {threads} \
					-i $PREFIX/asm.ctg.lay \
					-o {output}
		else
			touch {output}
		fi
	
	########################## MINIASM ##########################################
	# TODO: fix the output to generate a fasta
	elif [ {wildcards.ASM} == 'miniasm' ]; then 
		# Overlap for PacBio reads (or use "-x ava-ont" for nanopore read overlapping)
		minimap2 -x {overlaper} -t {threads} {input} {input} | \
				gzip -1 > $PREFIX/reads.paf.gz
		
		# Layout
		miniasm -f {input} $PREFIX/reads.paf.gz > $PREFIX/reads.gfa

		if [ -s $PREFIX/reads.gfa ] ; then 
			{base}/gfa2fasta.sh $PREFIX/reads.gfa {output}
		else
			touch {output}
		fi

	fi 

else
	>&2 echo " no real assembly"
	> {output}
fi
'''

# check if the assembly is not size 0
rule assemblyReport:
    input:  
        oasm= 'group.{n}/{ASM}.assembly/asm.contigs.fasta',
    output: 
        asm=  'group.{n}/{ASM}.assembly.fasta',
        report='group.{n}/{ASM}.report.txt'
    shell:
        """
        if [ -s {input.oasm} ]
        then
            cp {input.oasm} {output.asm}
	        echo "Number of reads " > {output.report}
	        echo "Assembly number of contigs" >> {output.report}
        else
            touch {output.asm}
            touch {output.report}
        fi
        """ 

if( ISONT or MM2 ):
	rule makeFakeCorrected:
		input:
			asm= 'group.{n}/{ASM}.assembly.fasta'
		output:
			quiver= 'group.{n}/{ASM}.assembly.consensus.fasta',
		threads: 4
		shell:
			"""
			cp {input.asm} {output.quiver}
			"""
else:
	rule bamFromAssembly:
		input:
			asm= 'group.{n}/{ASM}.assembly.fasta',
			H2= 'group.{n}/H2.WH.bam',
		output: 
			asmbam= 'group.{n}/{ASM}.assembly.bam',
		threads: 4
		shell:
			"""
			if [ ! -s {input.asm} ] 
			then
				# create empty files, this will allow other rules to conitnue 
				> {output.asmbam}
			else 
				blasr {input.H2} {input.asm} \
					--clipping subread --bam --bestn 1 --nproc {threads} \
					--out - | \
					samtools sort -m 4G -T tmp -o {output.asmbam}
				
				samtools index {output.asmbam}
			
			fi
			"""

	rule quiverFromBam:
		input:
			asmbam= 'group.{n}/{ASM}.assembly.bam',
			asm= 'group.{n}/{ASM}.assembly.fasta'
		output:
			quiver= 'group.{n}/{ASM}.assembly.consensus.fasta',
		threads: 4
		shell:
			'''
			# check 
			if [ ! -s {input.asm} ] 
			then
				# create empty files, this will allow other rules to conitnue 
				> {output.quiver}
			else
				source deactivate  2> /dev/null
				source {python2}
				which quiver
				quiver \
					--noEvidenceConsensusCall nocall --minCoverage 10 -j {threads} \
					-r {input.asm} -o {output.quiver} {input.asmbam} \
					|| ( >&2 echo " QUIVER FAILED TO RUN " && \
					cp {input.asm} {output.quiver} )
				
				# add the head of the non quivered file
				header=$(head -n 1 {input.asm})
				header2=">group.{wildcards.n}_quiver "$(echo $header | sed -e 's/>//')
				sed -i "1s/.*/$header2/" {output.quiver} 
			fi
			'''
#-----------------------------------------------------------------------------------------------------#






#--------------------------------------------------------------------------------------------#
# combines the output assemblies 
rule combineAsm:
	input:
		quiver= expand('group.{ID}/{{ASM}}.assembly.consensus.fasta', ID=IDS),
	output: 
		asm='{ASM}.assemblies.pre.pilon.fasta',
	run:
		collapse = os.path.basename(os.getcwd())
		rtn = ""
		counter = 1
		toAdd = []
		asmfile="{}.assembly.consensus.fasta".format(wildcards["ASM"])
		for asm in sorted( glob.glob("group.*/" + asmfile )): 
			match = re.match("(group.\d+)/" + asmfile, asm)
			group = match.group(1)
			print(group)
			recs = list(SeqIO.parse(asm, "fasta"))
			for rec in recs:
				rec.id = "{}_collapse.{}_id.{}".format(group, collapse, counter)
				rec.name = rec.id
				rec.seq = rec.seq.strip("N")
				counter += 1
				print(rec.id)
				toAdd.append(rec)
		#print(rtn)
		SeqIO.write(toAdd ,output["asm"], "fasta" ) 
	



# this needs more work to actually imporve thigns... right now it fixes nothing,
# need to use something not pilon. 
if(os.path.exists("illumina.orig.bam")):
	rule fastqIllumina:
		input:
			bam = "illumina.orig.bam",
		output:
			fastq = "illumina.fastq",
		shell:
			"""
			samtools bam2fq {input.bam} > {output.fastq}
			"""

	rule indexForBWA:
		input:
			asm='{ASM}.assemblies.pre.pilon.fasta',
		output:
			amb = "bwa/{ASM}.bwa.amb",
			ann = "bwa/{ASM}.bwa.ann",
			bwt = "bwa/{ASM}.bwa.bwt",
			pac = "bwa/{ASM}.bwa.pac",
			sa = "bwa/{ASM}.bwa.sa",
		shell:
			"""
			mkdir -p bwa
			bwa index -p bwa/{wildcards.ASM}.bwa -a bwtsw {input.asm}
			"""

	rule reMapIllumina:
		input:
			fastq = "illumina.fastq",
			asmWH='{ASM}.assemblies.pre.pilon.fasta',
			amb = "bwa/{ASM}.bwa.amb",
			ann = "bwa/{ASM}.bwa.ann",
			bwt = "bwa/{ASM}.bwa.bwt",
			pac = "bwa/{ASM}.bwa.pac",
			sa = "bwa/{ASM}.bwa.sa",

		output:
			bam = "{ASM}.illumina.bam",
		threads: 8
		shell:
			"""
			bwa mem -M -t {threads} -p bwa/{wildcards.ASM}.bwa {input.fastq} | \
					samtools view -bS - | \
					samtools sort - -o {output.bam}
			"""
			
	rule runPilon:
		input:
			asmWH='{ASM}.assemblies.pre.pilon.fasta',
			bam = "{ASM}.illumina.bam",
		output:
			asmWH = "{ASM}.assemblies.fasta",
			bai = "{ASM}.illumina.bam.bai",
		threads: 8
		shell:
			"""
			samtools index {input.bam}
			mkdir -p pilon_out
			pilon \
					--threads {threads} \
					--genome {input.asmWH} \
					--bam {input.bam} \
					--outdir pilon_out \
					--fix "indels" \
					--changes --vcf --tracks \
					--duplicates 
			cp pilon_out/pilon.fasta {output.asmWH}
			"""

else:
	rule copyToEnd:
		input:
			asm='{ASM}.assemblies.pre.pilon.fasta',
		output:
			asm = "{ASM}.assemblies.fasta",
		shell:
			"""
			cp {input.asm} {output.asm}
			"""


#
#-----------------------------------------------------------------------------------------------------#
#
if(os.path.exists("duplications.fasta")):
	rule truePSVs:
		input:
			dup="duplications.fasta",
			ref='ref.fasta',
			cuts='CC/mi.gml.cuts',
			vcf='snvs/assembly.consensus.nucfreq.vcf'
		output:
			refsam="truth/refVSdup.sam",
			refsnv="truth/refVSdup.snv",
			truthmatrix="truth/truth.matrix",
			sv="truth/refVSdup.SV"
		shell:
			"""
			mkdir -p truth
			
			#b {input.dup} {input.ref} --sam --bestn 1 --clipping subread > {output.refsam} 
			minimap2 \
				{input.ref} {input.dup} \
				-k 11 -a --eqx \
				-r {bandwidth} \
				-x asm20 | \
				samtools view -h -F 2308 - | \
				samtools sort -m 4G -T tmp -o {output.refsam}


			source {python2}

			{base}/PrintGaps.py \
				{input.ref} {output.refsam} --snv {output.refsnv} > {output.sv}

			{base}/CompareRefSNVWithPSV.py \
				--ref {output.refsnv} --psv {input.cuts} --vcf {input.vcf} \
				--refFasta {input.ref} --writevcf truth/true > {output.truthmatrix}

			"""

	rule truePSVsWithRefCordinates:
		input:
			refsnv="truth/refVSdup.snv",
		output:
			truth="truth/README.txt",
			snv = "truth/all_true.snv",
		run:
			# reads snv file into a dictoriry based on positon
			snvfile = open(input["refsnv"])
			allsnvs = {}
			for snvline in snvfile:
				token = snvline.split("\t")
				key = "{}_{}".format(token[0], token[2])
				allsnvs[key] = snvline.strip()	

			# reads all the truth files
			truesnvs=""
			for f in sorted(glob.glob("truth/true.*.vcf")):
				vcf = open(f)
				for line in vcf:
					token = line.split("\t")
					if(line[0]=="#" or len(token) < 2 ):
						continue
					key = "{}_{}".format(token[0], token[1])
					toadd = "{}\t{}\t{}\n".format(allsnvs[key], token[2], f)
					truesnvs += toadd
				vcf.close()
			open(output["snv"], "w+").write(truesnvs)

			shell('echo "exists" > {output.truth}')




	rule mapToRefAndDupsBlasr:
		input:	
			asm="{ASM}.assemblies.fasta",
			ref="ref.fasta",
			dup="duplications.fasta",
		output:
			refsam="asms/{ASM}.sam",
			dupsam="asms/{ASM}.dup.sam",
		threads: 8
		shell:"""
#b --nproc {threads} --sam --out - \
#	--bestn 1 --minMatch 11 --maxMatch 15 --nCandidates 50 \
#	{input.asm} {input.ref} | \
#	samtools view -h -F 4 - | samtools sort -m 4G -T tmp -o {output.refsam}

minimap2 \
	{input.ref} {input.asm} \
	-k 11 -a --eqx -t {threads} \
	-r {bandwidth} \
	-x asm20 | \
	samtools view -h -F 2308 - | \
	samtools sort -m 4G -T tmp -o {output.refsam}


#b --nproc {threads} --sam --out - \
#	--bestn 1 --minMatch 11 --maxMatch 15 --nCandidates 50 \
#	{input.asm} {input.dup} | \
#	samtools view -h -F 4 - | samtools sort -m 4G -T tmp -o {output.dupsam}

minimap2 \
	{input.dup} {input.asm} \
	-k 11 -a --eqx -t {threads} \
	-r {bandwidth} \
	-x asm20 | \
	samtools view -h -F 2308 - | \
	samtools sort -m 4G -T tmp -o {output.dupsam}
"""


	rule getTablesFromSam:
		input:
			refsam="asms/{ASM}.sam",
			dupsam="asms/{ASM}.dup.sam",
		output:
			dup="asms/{ASM}.dup.tbl",
			ref="asms/{ASM}.tbl",
		shell:
			"""
			{base}/samIdentity.py --header {input.refsam} > {output.ref}
			{base}/samIdentity.py --header {input.dupsam} > {output.dup}
			"""



	rule depthOnDuplications:
		input:
			reads="reads.sample.bam",
			ref = "duplications.fasta",
		output:
			bam = "asms/reads.dup.bam",
			depth="asms/dup_depth.tsv",
		threads:8
		shell:"""
#b {input.reads} {input.ref} \
#	--sam --nproc {threads} --out - \
#	--minAlignLength 500 \
#	--clipping subread | \
#	samtools view -bS -F 4 - | \
#	samtools sort -T tmp -o {output.bam}

samtools fasta {input.reads} |
	minimap2 \
		{input.ref} /dev/stdin \
		-k 11 -a --eqx -t {threads} \
		-r {bandwidth} \
		-x {MM2_TYPE} | \
			samtools view -bS -F 2308 - | \
			samtools sort -m 4G -T tmp -o {output.bam}

samtools index {output.bam}
samtools depth -aa {output.bam} > {output.depth}
"""


	rule plot_seqs_on_dup:
		input:
			depth="asms/dup_depth.tsv",
			sam="asms/{ASM}.dup.sam"
		output:
			plot = "{ASM}.SeqsOnDup.png",
		shell:
			"""
			{base}/plotDepth.py {input.depth} {output} --sam {input.sam}
			"""



	#
	# runs a summary script that just consilidates some data, which is later sued in plotting
	#
	rule summary:
		input:
			bed="ref.fasta.bed",
			plot = "{ASM}.SeqsOnDup.png",
			combine='{ASM}.assemblies.fasta',
			truth="truth/README.txt",
			tbldup=rules.getTablesFromSam.output.dup,
			tblref=rules.getTablesFromSam.output.ref,
			truthmatrix="truth/truth.matrix",
		output:
			summary="{ASM}.summary.txt",
			table="{ASM}.sda.table.tsv",
		shell:
			"""
			{base}/summary.py --assembler {wildcards.ASM} --summary {output.summary}
			{base}/overlapOfReadsCheck.py "group.*/H2.WH.bam" > read_collision.txt
			"""

		
	#
	#
	#
	rule bedForTrack:
		input:
			bedx="ref.fasta.bed",
			table="{ASM}.sda.table.tsv",
			summary="{ASM}.summary.txt",
			truthmatrix="truth/truth.matrix",
		output:
			asmbed="{ASM}.asm.bed",
		params:
			project=config["project"],
		shell:
			"""
			{base}/bedForABP.py --stats {input.table} \
					--summary {input.summary} --track {params.project} \
					--out {output.asmbed}
			"""

else:
	rule noDuplicaitonsFasta:
		input:
			ref='ref.fasta',
		output:
			truth="truth/README.txt",
			asmbed=expand("{ASM}.asm.bed", ASM=assemblers),
		shell:
			"""
			mkdir -p truth
			echo "does not exist" > {output.truth}
			touch {output.asmbed}
			"""
#-----------------------------------------------------------------------------------------------------#






#
# create a map of the coverage across the assembled duplications
#
rule coverageOnAsms:
	input:
		asm = "{ASM}.assemblies.fasta",
		reads = "reads.sample.bam",
		ref = "ref.fasta",
	output:
		refplus="asms/{ASM}.plus.Ref.fasta",
		cov="asms/{ASM}.depth.tsv",
		bam="asms/{ASM}.reads_on_asm.bam",
	threads:8
	shell:"""
cat {input.ref} {input.asm} > {output.refplus}

#b {input.reads} {output.refplus} --clipping subread \
#	--nproc {threads} --bestn 1 --sam --out - | \
#	samtools view -bS - | \
#	samtools sort -m 4G -o {output.bam} - 

samtools fasta {input.reads} |
	minimap2 \
		{output.refplus} /dev/stdin \
		-k 11 -a --eqx -t {threads} \
		-r {bandwidth} \
		-x {MM2_TYPE} | \
			samtools view -bS -F 2308 - | \
			samtools sort -m 4G -T tmp -o {output.bam}

samtools index {output.bam}
samtools depth -aa {output.bam} > {output.cov}
"""


rule plotCovOnAsm:
	input:
		cov="asms/{ASM}.depth.tsv",
	output:
		plot="{ASM}.CoverageOnAsm.png",
	shell:
		"""
		{base}/plotDepth.py {input.cov} {output.plot}
		"""




#-----------------------------------------------------------------------------------------------------#
if(os.path.exists("real.fasta")):
	# create a map of the reads onto the real end results 
	rule coverageOnReal:
		input:
			ref = "real.fasta",
			reads = "reads.fofn",
		output:
			cov="real/real_depth.tsv",
			bam="real/reads_on_real.bam",
		threads:16
		shell:
			"""
			blasr {input.reads} {input.ref} --bestn 1 --clipping subread --nproc {threads} --sam --out - | \
					samtools view -bS -F 4 - | \
					samtools sort -m 4G -o {output.bam} - 
			samtools index {output.bam}
			samtools depth -aa {output.bam} > {output.cov}
			"""
	rule map_asms_to_real:
		input:
			asm="{ASM}.assemblies.fasta",
			ref="real.fasta"
		output:
			m5="real/{ASM}.real.m5",
			sam="real/{ASM}.real.sam",
		shell:
			"""
			blasr -m 5 --bestn 1 --out {output.m5} {input.asm} {input.ref}
			blasr -m 5 --bestn 1 --sam --clipping subread --out {output.sam} {input.asm} {input.ref}
			"""


	rule plotCovOnReal:
		input:
			cov="real/real_depth.tsv",
			sam="real/real.sam",
		output:
			plot="real/SeqsOnReal.png",
			done="real/done.txt",
		shell:
			"""
			{utils}/plotDepth.py {input.cov} {output.plot} --sam {input.sam}
			touch {output.done}
			"""

else:
	rule FakeCovOnReal:
		input:
			asm="{ASM}.assemblies.fasta",
		output:
			done="real/{ASM}.done.txt",
		shell:
			"""
			touch {output.done}
			"""
#-----------------------------------------------------------------------------------------------------#




rule final:
	input:
		combine=expand('{ASM}.assemblies.fasta', ASM=assemblers ),
		plot=expand("{ASM}.CoverageOnAsm.png", ASM=assemblers),
		asmbed=expand("{ASM}.asm.bed", ASM=assemblers),
		done=expand("real/{ASM}.done.txt", ASM=assemblers),
		asms=expand("group.{ID}/{ASM}.assembly/asm.contigs.fasta", ID=IDS, ASM=assemblers),
		truth="truth/README.txt",
		dupindex="read.duplication.index",
	output: 'final'
	shell:"""
touch {output}
# remove extra files from assemblies, this speeds up the dag building for snakemake by a lot
rm -rf \
group.*/canu.assembly/canu-logs \
group.*/canu.assembly/canu-scripts \
group.*/canu.assembly/correction \
group.*/canu.assembly/correction.html.files \
group.*/canu.assembly/trimming \
group.*/canu.assembly/trimming.html.files \
group.*/canu.assembly/unitigging \
group.*/canu.assembly/unitigging.html.files 
"""

