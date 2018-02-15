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

blasr= '~mchaisso/projects/AssemblyByPhasing/scripts/abp/bin/blasr'
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
	"coverage.json"	
	
groups= glob.glob("group.[0-9]*.vcf")
IDS= []
for group in groups:
    group= group.strip().split(".")
    assert group[0]=="group"
    IDS.append(group[1])


rule master:
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
    input: "reads.bam"
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
    threads: 1 # for some reason loading this env multiple times breaks it, so stop from parallel exe
                # this should no longer be nessisary as whatshap is installed in the new anaconda
    shell:
        """
        source {base2}/env_whatshap.cfg 
        #source activate whatshap  # this should no longer be nessisary as whatshap is installed in the new anaconda
        whatshap haplotag -o  {output.hap} --reference ref.fasta {input.hapvcf} {input.hapbam} 
        samtools view -h -o - {output.hap} | grep -E "^@|HP:i:1" >  {output.hapH1} 
        samtools view -h -o - {output.hap} | grep -E "^@|HP:i:2" >  {output.hapH2}
        """
#-----------------------------------------------------------------------------------------------------#



#-----------------------------------------------------------------------------------------------------#
#
# This is also a partitioning script but it was written by mark, and does not seem to work quite as well
# however I am still running it in order to comapre results
#
rule partitionReads:
    input: 
        vcf=  'group.{n}.vcf',
        pgroup= 'group.{n}/group.{n}',
        bam=  'reads.bam',
        # the following just makes sure whatshap is run first
        hapH2= 'group.{n}/H2.WH.sam'
    output:
        # these are switched from 1 and 2 because marks partition is the opposite of whatshap and I want
        # them to be the same H1 H2 convention from here on out
        group1='group.{n}/H2.Mark.sam',
        group2='group.{n}/H1.Mark.sam'
    shell:
    	"""
        #samtools view -h {input.bam} > temp.txt
        samtools view -h {input.bam} \
	    	| ~mchaisso/projects/pbgreedyphase/partitionByPhasedSNVs \
		    	--vcf {input.vcf} \
			    --ref ref.fasta \
			    --h1 {output.group1}  --h2 {output.group2} --sam /dev/stdin \
			    --unassigned /dev/null \
				    	--phaseStats group.{wildcards.n}/group.stats \
					    --block 4 \
					    --minGenotyped 2 \
					    --minDifference 3 \
                        || touch {output.group1} && touch {output.group2}
        """

# generate a fasta file from the partition
rule fastaFromPartition:
    input: 
        group1='group.{n}/H1.Mark.sam',
        group2='group.{n}/H2.Mark.sam'
    output:
        pfasta1='group.{n}/H1.Mark.fasta',
        pfasta2='group.{n}/H2.Mark.fasta'
    shell:
        '''
        grep -v "^@" {input.group1}  | awk '{{ print ">"$1; print $10;}}' > {output.pfasta1} 
        grep -v "^@" {input.group2}  | awk '{{ print ">"$1; print $10;}}' > {output.pfasta2}
        '''

        #'make -f {base2}/PartitionReads.mak VCF={input.vcf} OUTDIR=group.{wildcards.n}/group HAP=2'
#-----------------------------------------------------------------------------------------------------#



#-----------------------------------------------------------------------------------------------------#
#
# This part runs the assembly based on the partitions 
#

# get fasta files from the reads
# this should be changed if we decide to drop group.1.sam
rule readsFromSam:
    input: 
        H2= 'group.{n}/H2.{prefix}.sam',
        # the following just makes sure whatshap was run and that marks partition was run
        #hapH2=  'group.{n}/H2.WH.sam',
        #markH2= 'group.{n}/H2.Mark.sam' 
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

		#cat group.{wildcards.n}/{wildcards.prefix}.temp.txt | {PBS}/local_assembly/StreamSamToFasta.py | \
        #    ~mchaisso/projects/PacBioSequencing/scripts/falcon/FormatFasta.py --fakename  > \
		#	{output.pfasta} || \
		#	touch {output.pfasta}


        rm -f group.{wildcards.n}/{wildcards.prefix}.temp.txt 
        """

# run the assembly
rule runAssembly:
    input: 'group.{n}/{prefix}.reads.fasta'
    output: 'group.{n}/{prefix}.assembly/asm.contigs.fasta'
    threads: 16
    shell:
        """
		# make sure we actaully re run the assembly
		rm -rf group.{wildcards.n}/{wildcards.prefix}.assembly/*

        if [ -s {input} ]; then
            module load java/8u25 && {CANU_DIR}/canu -pacbio-raw {input} genomeSize=60000 \
                -d group.{wildcards.n}/{wildcards.prefix}.assembly \
		        -p asm useGrid=false  gnuplotTested=true  corMhapSensitivity=high corMinCoverage=1 \
		        maxThreads={threads} cnsThreads={threads} ovlThreads={threads} \
		        mhapThreads={threads} contigFilter="2 1000 1.0 1.0 2" \
                || ( >&2 echo " no real assembly" && \
                mkdir -p group.{wildcards.n}/{wildcards.prefix}.assembly && \
                > {output} )

        else
            >&2 echo " no real assembly"
            mkdir -p group.{wildcards.n}/{wildcards.prefix}.assembly
            > {output}
        fi
        """

# this is currently not being run, uncooment the input line from assemblyReport, to start running it, and add additional code 
# to handle it
rule runFalcon:
    input:
        canu = 'group.{n}/{prefix}.assembly/asm.contigs.fasta',
        reads = 'group.{n}/{prefix}.reads.fasta',
    output:
       'group.{n}/{prefix}.assembly/falcon.readme'
    threads: 16
    shell:
        """
        # if the assembly is empty lets try out falcon, and the reads file is not empty 
        if [ ! -s  {input.canu} ] && [ -s {input.reads} ]; then
            
            make a falcon dir
            mkdir -p group.{wildcards.n}/{wildcards.prefix}.assembly/falcon

            # put the reads in an fofn for falcon
            echo $(readlink -f {input.reads}) > group.{wildcards.n}/{wildcards.prefix}.assembly/falcon/input.fofn
            
            # move into the assembly directory 
            pushd group.{wildcards.n}/{wildcards.prefix}.assembly/falcon/

            # setup falcon
            PBS=~mchaisso/projects/PacBioSequencing/scripts
            BLASR=~mchaisso/projects/blasr-repo/blasr
            MMAP=~mchaisso/software/minimap
            MASM=~mchaisso/software/miniasm
            QUIVER=~mchaisso/software/quiver
            PBG=~mchaisso/projects/pbgreedyphase
            
            # run falcon
            source ~mchaisso/scripts/setup_falcon.sh && fc_run.py ~mchaisso/projects/PacBioSequencing/scripts/local_assembly/falcon/fc_run.low_cov.cfg
           
            # move the assembly into the spot of the other assembly 
            cp 2-asm-falcon/p_ctg.fa ../asm.contigs.fasta

            popd 

            echo "Falcon Ran" > {output} 

        else
            # not running falcon
            echo "Canu worked so falcon did not run" > {output}
        fi
        
        """

# check if the assembly is not size 0
rule assemblyReport:
    input:  
        oasm= 'group.{n}/{prefix}.assembly/asm.contigs.fasta',
        preads='group.{n}/{prefix}.reads.fasta',
        #falcon='group.{n}/{prefix}.assembly/falcon.readme'
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
    threads: 16
    shell:
        """
        if [ ! -s {input.asm} ] 
	    then
            # create empty files, this will allow other rules to conitnue 
		    > {output.asmbam}
            > {output.asmbas}
        else 
            ~mchaisso/projects/blasr-repo/blasr/pbihdfutils/bin/samtobas {input.H2} {output.asmbas} -defaultToP6
	        ~mchaisso/projects/blasr-repo/blasr/alignment/bin/blasr {output.asmbas} {input.asm} \
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
    threads: 16
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




#-----------------------------------------------------------------------------------------------------#
#
# map back to duplications.fasta to determine the region wiht the highest %ID and 
# the average %Id acrosss all the regions in the human reference
# (in the furture I may add soemthing that does this for non quivered files)
#
if(os.path.exists("duplications.fasta")):
	rule newName:
		input:
			dup="duplications.fasta",
		output:
			dup="duplications.fixed.fasta",
		run:
			fixed = ""
			idx = 1
			f = open(input["dup"])
			for line in f:
				line = line.strip()
				if(line[0]==">"):
					line += "\tcopy" + str(idx)
					idx += 1
				fixed += line + "\n"

			open(output["dup"], "w+").write(fixed)
		'''
		shell:
			"""
			awk 'BEGIN{{count=1}}{{if($0~/^>/){{print ">copy"count,$0;count++}}else{{print}}}}' \
				 < duplications.fasta | sed -e 's/ >chr/\tchr/g' > {output.dup}
			"""
		'''

	rule bestMappings:
		input:
			dup="duplications.fixed.fasta",
			quiver= 'group.{n}/{prefix}.assembly.consensus.fasta',
			preads= 'group.{n}/H2.{prefix}.sam',
		output:
			best= 'group.{n}/{prefix}.best.m4',
			average= 'group.{n}/{prefix}.average.m4',
			best5= 'group.{n}/{prefix}.best.m5',
			average5= 'group.{n}/{prefix}.average.m5',
		shell:
			"""
			if [ ! -s {input.quiver} ] 
			then
			# create empty files, this will allow other rules to conitnue 
			> {output.best}
			> {output.average}
			> {output.best5}
			> {output.average5}
			else
			{blasr} {input.dup} {input.quiver} -bestn 1 -header -m 4 > {output.average}
			{blasr} {input.quiver} {input.dup} -bestn 1 -header -m 4 > {output.best}
			{blasr} {input.dup} {input.quiver} -bestn 1 -m 5 > {output.average5}
			{blasr} {input.quiver} {input.dup} -bestn 1 -m 5 > {output.best5}
			fi
			""" 

	rule mapBackPartition:
		input:
			dup="duplications.fixed.fasta",
			quiver= 'group.{n}/{prefix}.assembly.consensus.fasta',
			#preads= 'group.{n}/{prefix}.reads.fasta',
			preads= 'group.{n}/H2.{prefix}.sam',
		output:
			psam= 'group.{n}/{prefix}.{n}.sam',
			pbam= 'group.{n}/{prefix}.{n}.bam',
		shell:
			"""
			if [ ! -s {input.preads} ] 
			then
			# create empty files, this will allow other rules to conitnue 
			> {output.psam}
			> {output.pbam}
			else
			{blasr} {input.preads} {input.dup} -bestn 1 -sam > {output.psam}
			samtools view -b {output.psam} \
					| samtools sort -m 4G -o {output.pbam}
			#samtools view -S -b -o temp.bam {output.psam}
			#samtools sort -m 4G -o {output.pbam} temp.bam
			samtools index {output.pbam}
			#rm temp.bam
			fi
			"""

	rule truePSVs:
		input:
			dup="duplications.fixed.fasta",
			ref='ref.fasta',
			cuts='mi.gml.cuts',
			vcf='assembly.consensus.nucfreq.vcf'
		output:
			truth="truth/README.txt",
			refsam="refVSdup.sam",
			refsnv="refVSdup.snv",
			truthmatrix="truth.matrix"
		shell:
			"""
			mkdir -p truth
			echo "exists" > {output.truth}

			blasr {input.dup} {input.ref} -sam -bestn 1 -clipping soft > {output.refsam} 

			/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/PrintGaps.py \
				{input.ref} {output.refsam} --snv {output.refsnv} > refVSdup.SV

			~mchaisso/projects/AssemblyByPhasing/scripts/utils/CompareRefSNVWithPSV.py \
				--ref {output.refsnv} --psv {input.cuts} --vcf {input.vcf} \
				--refFasta {input.ref} --writevcf truth/true > truth.matrix

			"""

	rule truePSVsWithRefCordinates:
		input:
			truth="truth/README.txt",
			refsnv="refVSdup.snv",
		output:
			snv = "all_true.snv",
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

else:
    rule noMapping:
        input:
            quiver= 'group.{n}/{prefix}.assembly.consensus.fasta',
            reads= 'group.{n}/H2.{prefix}.sam',
        output:
            best= 'group.{n}/{prefix}.best.m4',
            average= 'group.{n}/{prefix}.average.m4',
            best5= 'group.{n}/{prefix}.best.m5',
            average5= 'group.{n}/{prefix}.average.m5',
        shell:
            """
            > {output.average}
            > {output.best}
            > {output.best5}
            > {output.average5}
            """

    rule noMappingBest:
        input:
            quiver= 'group.{n}/{prefix}.assembly.consensus.fasta',
            reads= 'group.{n}/H2.{prefix}.sam',
        output:
            sam= 'group.{n}/{prefix}.best.sam',
            bam= 'group.{n}/{prefix}.best.bam',
            psam= 'group.{n}/{prefix}.{n}.sam',
            pbam= 'group.{n}/{prefix}.{n}.bam',
        shell:
            """
            > {output.bam}
            > {output.sam}
            > {output.pbam}
            > {output.psam}
            """
    
    rule noTruePSVs:
        input:
            ref='ref.fasta',
            cuts='mi.gml.cuts',
            vcf='assembly.consensus.nucfreq.vcf'
        output:
            truth="truth/README.txt",
        shell:
            """
            mkdir -p truth
            echo "does not exist" > {output.truth}
            """



#-----------------------------------------------------------------------------------------------------#




#--------------------------------------------------------------------------------------------# checks to see if the assembly exists and if not removes the empty file if it does not and all other 
# empty files
# creats a group output file
# then creats a empty file singinaling we are done
# removeing empty files is actually a bad idea with snakemake becuase it will try to re run canu next time, when we know it will jsut fail
#
rule removeEmptyAsm:
    input:
        # commenting the following line should remove Marks partitioning script from the required assembly
        #eMark=expand( 'group.{ID}/Mark.best.m4', ID=IDS),
        eWH=expand( 'group.{ID}/WH.best.m4',   ID=IDS),
        e5WH=expand( 'group.{ID}/WH.best.m5',   ID=IDS),
        esam=expand( 'group.{ID}/WH.{ID}.sam',   ID=IDS),
        ebam=expand( 'group.{ID}/WH.{ID}.bam',   ID=IDS),
    output: 'removeEmpty'
    shell:
        """
        # removes any empty assemblies we may have created along the way 
        #find group.*/ -maxdepth 1 -size  0  | xargs -n 1 rm -f 
        touch {output}
        """

#cat group.*/Mark.assembly.consensus.fasta > {output.asmMark} || > {output.asmMark}
rule combineAsm:
	input:
		remove='removeEmpty', 
	output: 
		asmWH='WH.assemblies.fasta',
		#asmMark='Mark.assemblies.fasta'
	run:
		collapse = os.path.basename(os.getcwd())
		shell("rm " + input["remove"] )
		#shell('cat group.*/WH.assembly.consensus.fasta > ' + output["asmWH"])# + ' || > ' + output["asmWH"])
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
# cross match command
#
CM = """export PATH=$PATH:/net/eichler/vol2/local/inhousebin
cross_match -tags -alignments -masklevel 0 -minmatch 50 -bandwidth 250 """

"""-gap_ext -1 -gap_init -10 -penalty -1""" 
"""#blasr -bestn 1 -sam -clipping soft -out {output.sam} {input.asmWH} {input.ref}"""

useBlasr=True

if(useBlasr):
	rule mapToRefAndDupsBlasr:
		input:	
			asmWH="WH.assemblies.fasta",
			ref="ref.fasta",
			dup="duplications.fixed.fasta",
		output:
			refsam="WH.sam",
			dupsam="WH_dup.sam",
		threads: 8
		shell:
			"""
			blasr -nproc {threads} -sam -clipping soft -out /dev/stdout \
				-bestn 1 -minMatch 11 -maxMatch 15 -nCandidates 50 \
				{input.asmWH} {input.ref} | \
				samtools view -h -F 4 - | samtools sort -m 4G -T tmp -o {output.refsam}
			
			blasr -nproc {threads} -sam -clipping soft -out /dev/stdout \
				-bestn 1 -minMatch 11 -maxMatch 15 -nCandidates 50 \
				{input.asmWH} {input.dup} | \
				samtools view -h -F 4 - | samtools sort -m 4G -T tmp -o {output.dupsam}
			
			"""
	rule getTablesFromSam:
		input:
			refsam="WH.sam",
			dupsam="WH_dup.sam",
		output:
			dup="WH_dup.tsv",
			ref="WH.tsv",
		shell:
			"""
			~mvollger/projects/utility/samIdentity.py --header {input.refsam} > {output.ref}
			~mvollger/projects/utility/samIdentity.py --header {input.dupsam} > {output.dup}
			"""


else:
	rule map_asms_to_ref:
		input:
			asmWH="WH.assemblies.fasta",
			ref="ref.fasta"
		output:
			cm="WH.cm",
		shell:
			"""
			{CM} {input.asmWH} {input.ref}  > {output.cm}
			"""

	#blasr -bestn 1 -sam -clipping soft {input.asmWH} {input.ref} -out {output.blasr}
	#{CM} {input.asmWH} {input.ref}  > {output.cm}
	rule map_asms_to_duplicaitons:
		input:
			asmWH="WH.assemblies.fasta",
			ref="duplications.fixed.fasta"
		output:
			#cm="WH_dup.cm",
			#blasr="WH_dup.blasr.sam",
			#{CM} {input.asmWH} {input.ref}  > {output.cm}
			blast="WH_dup.blast",
		threads: 8
		shell:
			"""
			blatdb="~mvollger/assemblies/hg38/ucsc.hg38.no_alts.fasta.2bit"		
			# for a blast alignmnet 
			source ~mvollger/projects/builds/anaconda/activateConda3.sh
			blastn -task megablast \
					-db ~mvollger/assemblies/hg38/ucsc.hg38.no_alts.fasta \
					-query {input.asmWH} \
					-out {output.blast} \
					-word_size 50 \
					-parse_deflines \
					-num_threads {threads} \
					-num_alignments 1 \
					-outfmt "7 score sstrand qseqid sseqid pident slen qlen length qcovhsp qcovs nident mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq"
					#-max_target_seqs 1 \
					#-max_hsps 1 \
				
			"""
	#
	rule map_asms_to_ref_convert:
		input:
			cm="WH.cm",
		output:
			sam="WH.sam",
			tsv="WH.tsv",
		shell:
			"""
			# the sam pipe excludes none primary alignmnets 
			~mvollger/projects/cmpSeq/cmToSam.py {input.cm} | samtools view -h -F 256 >  {output.sam}
			~mvollger/projects/utility/samIdentity.py --header {output.sam} > {output.tsv}
			"""

	rule map_asms_to_duplicaitons_convert:
		input:
			blast="WH_dup.blast",
			#cm="WH_dup.cm",
		output:
			sam="WH_dup.sam",
			tsv="WH_dup.tsv",
		shell:
			"""
			~mvollger/projects/cmpSeq/cmToSam.py --blast {input.blast} | \
					samtools view -h -F 256 >  {output.sam}
			~mvollger/projects/utility/samIdentity.py --header {output.sam} > {output.tsv}
			"""




rule plot_seqs_on_dup:
    input:
        depth="dup_depth.tsv",
        sam="WH_dup.sam"
    output:
        "seqsOnDup.png",
    shell:
        """
        {utils}/plotDepth.py {input.depth} {output} --sam {input.sam}
        """
'''
rule plot_seqs_on_cov:
    input:
        depth="depth.tsv",
        sam="WH.sam",
    output:
        "seqs.png",
    shell:
        """
        {utils}/plotDepth.py {input.depth} {output} --sam {input.sam}
        """
'''

#
# create a map of the coverage across the assembled duplications
#
rule coverageOnAsms:
	input:
		asm = "WH.assemblies.fasta",
		ref = "ref.fasta",
		reads = "reads.fasta",
	output:
		cov="asm_depth.tsv",
		bam="reads_on_asm.bam",
		refWH="refAndWH.fasta",
	threads:8
	shell:
		"""
		cat {input.ref} {input.asm} > refAndWH.fasta
		{blasr} {input.reads} refAndWH.fasta -clipping soft \
				-nproc {threads} -bestn 1 -sam -out /dev/stdout | \
				samtools view -bS - | \
				samtools sort -m 4G -o {output.bam} - 
				samtools index {output.bam}
		samtools depth -aa {output.bam} > {output.cov}
		"""

rule plotCovOnAsm:
	input:
		cov="asm_depth.tsv",
	output:
		plot="CoverageOnAsm.png",
	shell:
		"""
		{utils}/plotDepth.py {input.cov} {output.plot}
		"""

#
# runs a summary script that just consilidates some data, which is later sued in plotting
#
rule summary:
	input:
		bed="ref.fasta.bed",
		combine='WH.assemblies.fasta',
		truth="truth/README.txt",
		tsv="WH_dup.tsv",
	output:
		summary="summary.txt",
		table="abp.table.tsv",
	shell:
		"""
		{base2}/summary.py
		"""
#
#
#
rule bedForTrack:
	input:
		bedx="ref.fasta.bed",
		table="abp.table.tsv",
		summary="summary.txt",
	output:
		asmbed="asm.bed",
	params:
		project=config["project"],
	shell:
		"""
		bedForABP.py {input.table} {params.project}
		"""

#
#
#
rule miropeats:
	input:
		refWH="refAndWH.fasta",
	output:
		miro = "refWH.miro.pdf",
	shell:
		"""
		miropeats -s 500 -onlyinter {input.refWH} > contig_compare
		if [ -f threshold500 ]
		then
			mv threshold500 refWH.miro.ps
			ps2pdf refWH.miro.ps
		else
			touch {output.miro}
		fi
		
		"""
# pdf='mi.cuts.gml.pdf',
rule final:
	input:
		combine='WH.assemblies.fasta',
		miro = "refWH.miro.pdf",
		plot="CoverageOnAsm.png",
		summary="summary.txt",
		truth="truth/README.txt",
		#seqPNG="seqs.png",
		dupPND="seqsOnDup.png",
		asmbed="asm.bed",
		plotReal="SeqsOnReal.png",
		snv = "all_true.snv",
	output: 'PSV2_done'
	shell:
		"""
		touch {output}
		"""


#-----------------------------------------------------------------------------------------------------#
if(os.path.exists("real.fasta")):
	# create a map of the reads onto the real end results 
	rule coverageOnReal:
		input:
			ref = "real.fasta",
			reads = "reads.fofn",
		output:
			cov="real_depth.tsv",
			bam="reads_on_real.bam",
		threads:16
		shell:
			"""
			{blasr} {input.reads} {input.ref} -bestn 1 -clipping soft -nproc {threads} -sam -out /dev/stdout \
					| samtools view -bS -F 4 - | \
					samtools sort -m 4G -o {output.bam} - 
			
			#source /net/eichler/vol24/projects/sequencing/pacbio/software/smrtanalysis/current/etc/setup.sh \
			#		&& /net/eichler/vol24/projects/sequencing/pacbio/smrt-link/smrtcmds/bin/pbalign \
			#		{input.reads} {input.ref} {output.bam} \
			#		--nproc {threads} --algorithmOptions="--minRawSubreadScore 800 --bestn 1"
			
			samtools index {output.bam}
			samtools depth -aa {output.bam} > {output.cov}
			"""
	rule map_asms_to_real:
		input:
			asmWH="WH.assemblies.fasta",
			ref="real.fasta"
		output:
			WHm5="real.m5",
			sam="real.sam",
		shell:
			"""
			blasr -m 5 -bestn 1 -out {output.WHm5} {input.asmWH} {input.ref}
			blasr -m 5 -bestn 1 -sam -clipping soft -out {output.sam} {input.asmWH} {input.ref}
			"""


	rule plotCovOnReal:
		input:
			cov="real_depth.tsv",
			sam="real.sam",
		output:
			plot="SeqsOnReal.png",
		shell:
			"""
			{utils}/plotDepth.py {input.cov} {output.plot} --sam {input.sam}
			"""

else:
	rule FakeCovOnReal:
		input:
			asmWH="WH.assemblies.fasta",
		output:
			plot="SeqsOnReal.png",
		shell:
			"""
			touch {output.plot}
			"""


    
    
