import os
import glob
import re

SNAKEMAKE_DIR= os.path.dirname(workflow.snakefile)

shell.executable("/bin/bash")
#shell.prefix("source %s/env_PSV.cfg; set -eo pipefail; " % SNAKEMAKE_DIR)
shell.prefix("source %s/env_PSV.cfg; " % SNAKEMAKE_DIR)
#shell.suffix("2> /dev/null")

blasr= '~mchaisso/projects/AssemblyByPhasing/scripts/abp/bin/blasr'
blasrDir= '~mchaisso/projects/blasr-repo/blasr'
scriptsDir= '/net/eichler/vol5/home/mchaisso/projects/AssemblyByPhasing/scripts/abp'
#base2="/net/eichler/vol21/projects/bac_assembly/nobackups/scripts"
base2="/net/eichler/vol2/home/mvollger/projects/abp"
PBS="/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts" 
CADIR="/net/eichler/vol5/home/mchaisso/software/wgs-8.2/Linux-amd64/bin/"
CANU_DIR="/net/eichler/vol5/home/mchaisso/software/canu/Linux-amd64/bin"
CORES="4"

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
        'mkdir -p {output.group} && '
        'echo "{output.tag}" >  {output.tag} '
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
        'source {base2}/env_whatshap.cfg && '
        '{base2}/fixVCF.py --out {output} --vcf {input.vcf}'

# whatshap requires unique bam entries, make those
# also requires that pysam is installed, should be taken care of by loading the whatshap anaconda env 
rule bamForWhatsHap:
    input: "reads.bam"
    output: 'reads.sample.bam'
    shell:
        'source {base2}/env_whatshap.cfg && '
        '{base2}/changeBamName.py'

# index reads.sample.bam
rule indexWhatshapReads:
    input:'reads.sample.bam'
    output: 'reads.sample.bam.bai',
    shell:
        'samtools index {input} {output}'

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
        ' source {base2}/env_whatshap.cfg && '
        ' source activate whatshap && '
        ' whatshap haplotag -o  {output.hap} --reference ref.fasta {input.hapvcf} {input.hapbam} &&'
        ' samtools view -h -o - {output.hap} | grep -E "^@|HP:i:1" >  {output.hapH1} && '
        ' samtools view -h -o - {output.hap} | grep -E "^@|HP:i:2" >  {output.hapH2} '
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
    	'''
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
        '''
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
        '''
        grep -v "^@" {input.H2} > temp.txt
        if [ -s temp.txt ]
        then
            cat temp.txt | {PBS}/local_assembly/StreamSamToFasta.py | \
                ~mchaisso/projects/PacBioSequencing/scripts/falcon/FormatFasta.py --fakename  > {output.pfasta}
        else
            >&2 echo " no real fasta file for assembly"
            touch {output.pfasta}
        fi
        rm -f temp.txt
        '''

# run the assembly
rule runAssembly:
    input: 'group.{n}/{prefix}.reads.fasta'
    output: 'group.{n}/{prefix}.assembly/asm.contigs.fasta'
    threads: 16
    shell:
        '''
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
        '''

# check if the assembly is not size 0
rule assemblyReport:
    input:  
        oasm= 'group.{n}/{prefix}.assembly/asm.contigs.fasta',
        preads='group.{n}/{prefix}.reads.fasta'
    output: 
        asm=  'group.{n}/{prefix}.assembly.fasta',
        report='group.{n}/{prefix}.report.txt'
    shell:
        '''
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
        '''

rule bamFromAssembly:
    input:
        asm= 'group.{n}/{prefix}.assembly.fasta',
        H2= 'group.{n}/H2.{prefix}.sam'
    output: 
        asmbam= 'group.{n}/{prefix}.assembly.bam',
        asmbas= 'group.{n}/{prefix}.reads.bas.h5'
    threads: 16
    shell:
        '''
        if [ ! -s {input.asm} ] 
	    then
            # create empty files, this will allow other rules to conitnue 
		    > {output.asmbam}
            > {output.asmbas}
        else 
            ~mchaisso/projects/blasr-repo/blasr/pbihdfutils/bin/samtobas {input.H2} {output.asmbas} -defaultToP6
	        ~mchaisso/projects/blasr-repo/blasr/alignment/bin/blasr {output.asmbas} {input.asm} \
                -clipping subread -sam -bestn 1 -out /dev/stdout  -nproc {threads} \
                | samtools view -bS - | samtools sort -T tmp -o {output.asmbam}
	        samtools index {output.asmbam}
        fi
        '''

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
        output:
            dup="duplications.fixed.fasta",
        shell:
            '''
             awk 'BEGIN{{count=1}}{{if($0~/^>/){{print ">copy"count,$0;count++}}else{{print}}}}' \
                     < duplications.fasta | sed -e 's/ >chr/\tchr/g' > {output.dup}
            '''


    rule bestMappings:
        input:
            dup="duplications.fixed.fasta",
            quiver= 'group.{n}/{prefix}.assembly.consensus.fasta',
            preads= 'group.{n}/H2.{prefix}.sam',
        output:
            best= 'group.{n}/{prefix}.best.m4',
            average= 'group.{n}/{prefix}.average.m4',
        shell:
            '''
            if [ ! -s {input.quiver} ] 
	        then
                # create empty files, this will allow other rules to conitnue 
		        > {output.best}
                > {output.average}
            else
                {blasr} {input.dup} {input.quiver} -bestn 1 -header -m 4 > {output.average}
                {blasr} {input.quiver} {input.dup} -bestn 1 -header -m 4 > {output.best}
            fi
            '''

    rule mapBackPartition:
        input:
            dup="duplications.fixed.fasta",
            quiver= 'group.{n}/{prefix}.assembly.consensus.fasta',
            preads= 'group.{n}/H2.{prefix}.sam',
        output:
            psam= 'group.{n}/{prefix}.{n}.sam',
            pbam= 'group.{n}/{prefix}.{n}.bam',
        shell:
            '''
            if [ ! -s {input.preads} ] 
	        then
                # create empty files, this will allow other rules to conitnue 
                > {output.psam}
                > {output.pbam}
            else
                {blasr} {input.preads} {input.dup} -bestn 1 -sam > {output.psam}
                samtools view -S -b -o temp.bam {output.psam}
                samtools sort -o {output.pbam} temp.bam
                samtools index {output.pbam}
                rm temp.bam
            fi
            '''
    
    rule truePSVs:
        input:
            dup="duplications.fixed.fasta",
            ref='ref.fasta',
            cuts='mi.gml.cuts',
            vcf='assembly.consensus.nucfreq.vcf'
        output:
            truth="truth/README.txt",
            refsam="refVSdup.sam",
            refsnv="refVSdup.snv"
        shell:
            '''
            mkdir -p truth
            echo "exists" > {output.truth}
             
            blasr {input.dup} {input.ref} -sam -bestn 1 > {output.refsam} 
            
            /net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/PrintGaps.py \
                    {input.ref} {output.refsam} --snv {output.refsnv} > refVSdup.SV
           
            ~mchaisso/projects/AssemblyByPhasing/scripts/utils/CompareRefSNVWithPSV.py \
                    --ref {output.refsnv} --psv {input.cuts} --vcf {input.vcf} \
                    --refFasta {input.ref} --writevcf truth/true
            
            '''

else:
    rule noMapping:
        input:
            quiver= 'group.{n}/{prefix}.assembly.consensus.fasta',
            reads= 'group.{n}/H2.{prefix}.sam',
        output:
            best= 'group.{n}/{prefix}.best.m4',
            average= 'group.{n}/{prefix}.average.m4',
            sam= 'group.{n}/{prefix}.best.sam',
            bam= 'group.{n}/{prefix}.best.bam',
        shell:
            '''
            > {output.average}
            > {output.best}
            > {output.bam}
            > {output.sam}
            '''

    rule noMappingBest:
        input:
            quiver= 'group.{n}/{prefix}.assembly.consensus.fasta',
            reads= 'group.{n}/H2.{prefix}.sam',
        output:
            sam= 'group.{n}/{prefix}.best.sam',
            bam= 'group.{n}/{prefix}.best.bam',
        shell:
            '''
            > {output.bam}
            > {output.sam}
            '''
    
    rule noTruePSVs:
        input:
            ref='ref.fasta',
            cuts='mi.gml.cuts',
            vcf='assembly.consensus.nucfreq.vcf'
        output:
            truth="truth/README.txt",
        shell:
            '''
            mkdir -p truth
            echo "does not exist" > {output.truth}
            '''



#-----------------------------------------------------------------------------------------------------#




#-----------------------------------------------------------------------------------------------------#
#
# checks to see if the assembly exists and if not removes the empty file if it does not and all other 
# empty files
# creats a group output file
# then creats a empty file singinaling we are done
#
rule removeEmptyAsm:
    input:
        eMark=expand( 'group.{ID}/Mark.best.m4', ID=IDS),
        eWH=expand( 'group.{ID}/WH.best.m4',   ID=IDS),
        esam=expand( 'group.{ID}/WH.{ID}.sam',   ID=IDS),
        ebam=expand( 'group.{ID}/WH.{ID}.bam',   ID=IDS),
    output: 'removeEmpty'
    shell:
        # removes any empty assemblies we may have created along the way 
        'find group.*/ -maxdepth 1 -size  0  | xargs -n 1 rm -f      && '
        'touch {output}'

rule combineAsm:
    input:
        remove='removeEmpty', 
    output: 
        asmWH='WH.assemblies.fasta',
        asmMark='Mark.assemblies.fasta'
    shell:
        '''
        rm {input.remove}
        cat group.*/Mark.assembly.consensus.fasta > {output.asmMark} || > {output.asmMark}
        cat   group.*/WH.assembly.consensus.fasta > {output.asmWH}   || > {output.asmWH}
        '''
#
# runs a summary script that just consilidates some data, which is later sued in plotting
#
rule summary:
    input:
        combine='WH.assemblies.fasta',
    output:
        summary="summary.txt",
    shell:
        '''
        {base2}/summary.py
        '''

#
# makes a gephi version of the plot, 
# much nice imo 
#
rule gephi:
    input:
        combine='WH.assemblies.fasta'
    output:
        pdf='mi.cuts.gml.pdf'
    shell:
        '''
        {base2}/gephi/gephi.sh
        '''


rule final:
    input:
        combine='WH.assemblies.fasta',
        summary="summary.txt",
        pdf='mi.cuts.gml.pdf',
        truth="truth/README.txt",
    output: 'PSV2_done'
    shell:
        '''
        touch {output}
        '''
#-----------------------------------------------------------------------------------------------------#



    
    
