import os
import glob
import re
from Bio import SeqIO
import re

SNAKEMAKE_DIR= os.path.dirname(workflow.snakefile)

shell.executable("/bin/bash")
#shell.prefix("source %s/env_PSV.cfg; set -eo pipefail; " % SNAKEMAKE_DIR)
shell.prefix("source %s/env_PSV.cfg; set -eo pipefail" % SNAKEMAKE_DIR)
#shell.suffix("2> /dev/null")

#blasr= '~mchaisso/projects/AssemblyByPhasing/scripts/abp/bin/blasr'
blasrDir= '~mchaisso/projects/blasr-repo/blasr'
scriptsDir= '/net/eichler/vol5/home/mchaisso/projects/AssemblyByPhasing/scripts/abp'
#base2="/net/eichler/vol21/projects/bac_assembly/nobackups/scripts"
base2="/net/eichler/vol2/home/mvollger/projects/abp"
utils="/net/eichler/vol2/home/mvollger/projects/utility"
PBS="/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts" 
# Canu 1.5 seems to have muchhhh better performance over canu 1.6
CANU_DIR="/net/eichler/vol5/home/mchaisso/software/canu/Linux-amd64/bin"
#CANU_DIR="/net/eichler/vol21/projects/bac_assembly/nobackups/canu/Linux-amd64/bin"
#CANU_DIR="/net/eichler/vol2/home/mvollger/projects/builds/canu/Linux-amd64/bin"

configfile:
	"abp.config.json"	
blasr = config["blasr"]
# min cov is used to detemrine the filter for getting rid of low read assemblies. 
MINREADS = int(config["MINCOV"]*1.0/2.0)


groups= glob.glob("group.[0-9]*.vcf")
IDS= []
for group in groups:
    group= group.strip().split(".")
    assert group[0]=="group"
    IDS.append(group[1])


rule all:
    input: "PSV2_done"
    message: "Running PSV2"


# global wild card constraint on n whihc is the group idenitifier
wildcard_constraints:
    n="\d+"

#-----------------------------------------------------------------------------------------------------#
#
# make the group directories, and create a file that jsut acts as a tag for when the dir was created
#
rule makeGroupDirs:
    input:  expand('group.{ID}.vcf', ID=IDS)
    output: 
        group='group.{n}/',
        tag='group.{n}/group.{n}'
    shell:
        """
        rm -rf summary.py WH.assemblies.fasta
		rm -rf {output.group} # the assemblies will not re run properlly unless it starts fresh 
        mkdir -p {output.group} 
        echo "{output.tag}" >  {output.tag} 
        """
#-----------------------------------------------------------------------------------------------------#




#-----------------------------------------------------------------------------------------------------#
#
# this part partitions the reads based on the vcf files (PSVs)
#

# make a phased vcf for whatshap, just a format change
# also requires that pysam is installed, should be taken care of by loading the whatshap anaconda env 
rule phasedVCF:
    input: 'group.{n}/group.{n}',
        vcf= 'group.{n}.vcf'
    output: 'group.{n}/phased.{n}.vcf'
    shell:
        """
        source {base2}/env_whatshap.cfg 
        {base2}/fixVCF.py --out {output} --vcf {input.vcf}
        """

# whatshap requires unique bam entries, make those
# also requires that pysam is installed, should be taken care of by loading the whatshap anaconda env 
rule bamForWhatsHap:
    input: "reads.bam", expand('group.{ID}.vcf', ID=IDS)
    output: 'reads.sample.bam'
    shell:
        """
        source {base2}/env_whatshap.cfg 
        {base2}/changeBamName.py
        """

# index reads.sample.bam
rule indexWhatshapReads:
    input:'reads.sample.bam'
    output: 'reads.sample.bam.bai',
    shell:
        """
        echo {output}
        samtools index {input}
        """ 

# run whats hap and get the partitioned sam files
rule whatsHap:
    input: 
        hapbam= 'reads.sample.bam',
        hapbai= 'reads.sample.bam.bai',
        hapvcf= 'group.{n}/phased.{n}.vcf'
    output: 
        hap= 'group.{n}/haplotagged.bam',
        hapH1= 'group.{n}/H1.WH.sam',
        hapH2= 'group.{n}/H2.WH.sam'
    shell:
        """
        source {base2}/env_whatshap.cfg 
        #whatshap haplotag -o  {output.hap} --reference ref.fasta {input.hapvcf} {input.hapbam} 
        # the above is probably a bad idea, it causes a lot more reads to exist in more than one haplotype 
		whatshap haplotag -o  {output.hap} {input.hapvcf} {input.hapbam} 
        samtools view -h -o - {output.hap} | grep -E "^@|HP:i:1" >  {output.hapH1} 
        samtools view -h -o - {output.hap} | grep -E "^@|HP:i:2" >  {output.hapH2}
		
		
		# this check how much overlap there is of reads between groups  
		# overlapOfReadsCheck.py "group.*/H2.WH.sam" > read_collision.txt
        """
#-----------------------------------------------------------------------------------------------------#



#
#-----------------------------------------------------------------------------------------------------#
#
# This part runs the assembly based on the partitions 
#

# get fasta files from the reads
# this should be changed if we decide to drop group.1.sam
rule readsFromSam:
    input: 
        H2= 'group.{n}/H2.{prefix}.sam',
    output:
        pfasta= 'group.{n}/{prefix}.reads.fasta'
    shell:
        """
        grep -v "^@" {input.H2} > group.{wildcards.n}/{wildcards.prefix}.temp.txt \
			|| touch group.{wildcards.n}/{wildcards.prefix}.temp.txt

		if [ -s group.{wildcards.n}/{wildcards.prefix}.temp.txt ]; then
           cat group.{wildcards.n}/{wildcards.prefix}.temp.txt | {PBS}/local_assembly/StreamSamToFasta.py | \
               ~mchaisso/projects/PacBioSequencing/scripts/falcon/FormatFasta.py --fakename  > \
				{output.pfasta};
        else
            >&2 echo " no real assembly"
            touch {output.pfasta};
        fi

        rm -f group.{wildcards.n}/{wildcards.prefix}.temp.txt 
        """

# run the assembly
rule runAssembly:
    input: 'group.{n}/{prefix}.reads.fasta'
    output: 'group.{n}/{prefix}.assembly/asm.contigs.fasta'
    threads: 4
    shell:
        """
		# make sure we actaully re run the assembly
		rm -rf group.{wildcards.n}/{wildcards.prefix}.assembly/*

        if [ -s {input} ]; then
            module load java/8u25 && {CANU_DIR}/canu -pacbio-raw {input} \
				genomeSize=60000 \
				corOutCoverage=300 \
				corMhapSensitivity=high \
				corMinCoverage=1 \
				gnuplotTested=true  \
		        -p asm useGrid=false  \
                -d group.{wildcards.n}/{wildcards.prefix}.assembly \
		        maxThreads={threads} cnsThreads={threads} ovlThreads={threads} \
		        mhapThreads={threads} \
				contigFilter="{MINREADS} 5000 1.0 .75 {MINREADS}" \
                || ( >&2 echo " no real assembly" && \
                mkdir -p group.{wildcards.n}/{wildcards.prefix}.assembly && \
                > {output} )

        else
            >&2 echo " no real assembly"
            mkdir -p group.{wildcards.n}/{wildcards.prefix}.assembly
            > {output}
        fi
        """


# check if the assembly is not size 0
rule assemblyReport:
    input:  
        oasm= 'group.{n}/{prefix}.assembly/asm.contigs.fasta',
        preads='group.{n}/{prefix}.reads.fasta',
    output: 
        asm=  'group.{n}/{prefix}.assembly.fasta',
        report='group.{n}/{prefix}.report.txt'
    shell:
        """
        if [ -s {input.oasm} ]
        then
            cp {input.oasm} {output.asm}
	        echo "Number of reads " > {output.report}
	        grep -c ">" {input.preads} >> {output.report}
	        echo "Assembly number of contigs" >> {output.report}
	        module load numpy/latest; ~mchaisso/software/mcsrc/UTILS/pcl {input.preads} \
                    | ~mchaisso/scripts/stats.py >> {output.report}
	        rm -rf templocal
        else
            touch {output.asm}
            touch {output.report}
        fi
        """ 


# if samtobas fails it is becasue we started with no qual information and it kills it
rule bamFromAssembly:
    input:
        asm= 'group.{n}/{prefix}.assembly.fasta',
        H2= 'group.{n}/H2.{prefix}.sam', # samtobas is not robust enough to work if I start without qual info so using fasta
        preads='group.{n}/{prefix}.reads.fasta',
    output: 
        asmbam= 'group.{n}/{prefix}.assembly.bam',
        asmbas= 'group.{n}/{prefix}.reads.bas.h5'
    threads: 4
    shell:
        """
        if [ ! -s {input.asm} ] 
	    then
            # create empty files, this will allow other rules to conitnue 
		    > {output.asmbam}
            > {output.asmbas}
        else 
            ~mchaisso/projects/blasr-repo/blasr/pbihdfutils/bin/samtobas {input.H2} {output.asmbas} -defaultToP6
			{blasr} {output.asmbas} {input.asm} \
                -clipping subread -sam -bestn 1 -out /dev/stdout  -nproc {threads} \
                | samtools view -bS - | samtools sort -m 4G -T tmp -o {output.asmbam}
	        
            samtools index {output.asmbam}
        fi
        """

rule quiverFromBam:
    input:
        asmbam= 'group.{n}/{prefix}.assembly.bam',
        asm= 'group.{n}/{prefix}.assembly.fasta'
    output:
        quiver= 'group.{n}/{prefix}.assembly.consensus.fasta',
    threads: 4
    shell:
        '''
        # check 
        if [ ! -s {input.asm} ] 
	    then
            # create empty files, this will allow other rules to conitnue 
		    > {output.quiver}
        else
            samtools faidx {input.asm}
	        source ~mchaisso/software/quiver/setup_quiver.sh; ~mchaisso/software/quiver/bin/quiver \
		        --noEvidenceConsensusCall nocall --minCoverage 10 -j {threads} \
		        -r {input.asm} -o {output.quiver} {input.asmbam}
            
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
		quiver= expand('group.{ID}/WH.assembly.consensus.fasta', ID=IDS),
	output: 
		asmWH='WH.assemblies.fasta',
	run:
		collapse = os.path.basename(os.getcwd())
		rtn = ""
		counter = 1
		toAdd = []
		for asm in sorted( glob.glob("group.*/WH.assembly.consensus.fasta")): 
			match = re.match( "(group.\d+)/WH.assembly.consensus.fasta", asm)
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
		SeqIO.write(toAdd ,output["asmWH"], "fasta" ) 



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

			blasr {input.dup} {input.ref} -sam -bestn 1 -clipping subread > {output.refsam} 

			/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/PrintGaps.py \
				{input.ref} {output.refsam} --snv {output.refsnv} > {output.sv}

			~mchaisso/projects/AssemblyByPhasing/scripts/utils/CompareRefSNVWithPSV.py \
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
			asmWH="WH.assemblies.fasta",
			ref="ref.fasta",
			dup="duplications.fasta",
		output:
			refsam="asms/WH.sam",
			dupsam="asms/WH_dup.sam",
		threads: 8
		shell:
			"""
			{blasr} -nproc {threads} -sam -clipping subread -out /dev/stdout \
				-bestn 1 -minMatch 11 -maxMatch 15 -nCandidates 50 \
				{input.asmWH} {input.ref} | \
				samtools view -h -F 4 - | samtools sort -m 4G -T tmp -o {output.refsam}
			
			{blasr} -nproc {threads} -sam -clipping subread -out /dev/stdout \
				-bestn 1 -minMatch 11 -maxMatch 15 -nCandidates 50 \
				{input.asmWH} {input.dup} | \
				samtools view -h -F 4 - | samtools sort -m 4G -T tmp -o {output.dupsam}
			
			"""

	rule getTablesFromSam:
		input:
			refsam="asms/WH.sam",
			dupsam="asms/WH_dup.sam",
		output:
			dup="asms/WH_dup.tsv",
			ref="asms/WH.tsv",
		shell:
			"""
			~mvollger/projects/utility/samIdentity.py --header {input.refsam} > {output.ref}
			~mvollger/projects/utility/samIdentity.py --header {input.dupsam} > {output.dup}
			"""



	rule depthOnDuplications:
		input:
			reads="reads.fasta",
			ref = "duplications.fasta",
		output:
			bam = "asms/reads.dup.bam",
			depth="asms/dup_depth.tsv",
		threads:8
		shell:
			"""
			{blasr} -sam  \
					-nproc {threads} -out /dev/stdout \
					-minAlignLength 500 -preserveReadTitle -clipping subread \
					{input.reads} {input.ref} | \
					samtools view -bSh -F 4 - | \
					samtools sort -T tmp -o {output.bam}
			
			samtools depth -aa {output.bam} > {output.depth}
			"""


	rule plot_seqs_on_dup:
		input:
			depth="asms/dup_depth.tsv",
			sam="asms/WH_dup.sam"
		output:
			plot = "SeqsOnDup.png",
		shell:
			"""
			module purge
			. /etc/profile.d/modules.sh
			module load modules modules-init modules-gs/prod modules-eichler
			module load anaconda/20161130
			{utils}/plotDepth.py {input.depth} {output} --sam {input.sam}
			"""



	#
	# runs a summary script that just consilidates some data, which is later sued in plotting
	#
	rule summary:
		input:
			bed="ref.fasta.bed",
			plot = "SeqsOnDup.png",
			combine='WH.assemblies.fasta',
			truth="truth/README.txt",
			tsv="asms/WH_dup.tsv",
			truthmatrix="truth/truth.matrix",
		output:
			summary="summary.txt",
			table="abp.table.tsv",
		shell:
			"""
			{base2}/summary.py
			overlapOfReadsCheck.py "group.*/H2.WH.sam" > read_collision.txt
			"""

		
	#
	#
	#
	rule bedForTrack:
		input:
			bedx="ref.fasta.bed",
			table="abp.table.tsv",
			summary="summary.txt",
			truthmatrix="truth/truth.matrix",
		output:
			asmbed="asm.bed",
		params:
			project=config["project"],
		shell:
			"""
			bedForABP.py {input.table} {params.project}
			"""

else:
	rule noDuplicaitonsFasta:
		input:
			ref='ref.fasta',
		output:
			truth="truth/README.txt",
			asmbed="asm.bed",
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
		asm = "WH.assemblies.fasta",
		ref = "ref.fasta",
		reads = "reads.fasta",
	output:
		cov="asms/asm_depth.tsv",
		bam="asms/reads_on_asm.bam",
		refWH="asms/refAndWH.fasta",
	threads:8
	shell:
		"""
		cat {input.ref} {input.asm} > {output.refWH}
		{blasr} {input.reads} {output.refWH} -clipping subread \
				-nproc {threads} -bestn 1 -sam -out /dev/stdout | \
				samtools view -bS - | \
				samtools sort -m 4G -o {output.bam} - 
				samtools index {output.bam}
		samtools depth -aa {output.bam} > {output.cov}
		"""

rule plotCovOnAsm:
	input:
		cov="asms/asm_depth.tsv",
	output:
		plot="CoverageOnAsm.png",
	shell:
		"""
		module purge
		. /etc/profile.d/modules.sh
		module load modules modules-init modules-gs/prod modules-eichler
		module load anaconda/20161130
		{utils}/plotDepth.py {input.cov} {output.plot}
		"""

#
#
#
rule miropeats:
	input:
		refWH="asms/refAndWH.fasta",
	output:
		miro = "asms/refWH.miro.pdf",
	shell:
		"""
		miropeats -s 500 -onlyinter {input.refWH} 
		if [ -f threshold500 ]
		then
			mv threshold500 temp.ps
			ps2pdf temp.ps {output.miro}
			rm temp.ps
		else
			touch {output.miro}
		fi
		
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
			{blasr} {input.reads} {input.ref} -bestn 1 -clipping subread -nproc {threads} -sam -out /dev/stdout \
					| samtools view -bS -F 4 - | \
					samtools sort -m 4G -o {output.bam} - 
			samtools index {output.bam}
			samtools depth -aa {output.bam} > {output.cov}
			"""
	rule map_asms_to_real:
		input:
			asmWH="WH.assemblies.fasta",
			ref="real.fasta"
		output:
			WHm5="real/real.m5",
			sam="real/real.sam",
		shell:
			"""
			{blasr} -m 5 -bestn 1 -out {output.WHm5} {input.asmWH} {input.ref}
			{blasr} -m 5 -bestn 1 -sam -clipping subread -out {output.sam} {input.asmWH} {input.ref}
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
			asmWH="WH.assemblies.fasta",
		output:
			done="real/done.txt",
		shell:
			"""
			touch {output.done}
			"""
#-----------------------------------------------------------------------------------------------------#





rule final:
	input:
		combine='WH.assemblies.fasta',
		miro = "asms/refWH.miro.pdf",
		plot="CoverageOnAsm.png",
		truth="truth/README.txt",
		asmbed="asm.bed",
		done="real/done.txt",
	output: 'PSV2_done'
	shell:
		"""
		touch {output}
		"""



    
    
