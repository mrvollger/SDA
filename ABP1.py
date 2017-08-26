import os
import glob
from Bio import SeqIO
import re

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)

shell.executable("/bin/bash")
#shell.prefix("source %s/env_PSV.cfg; set -eo pipefail; " % SNAKEMAKE_DIR)
shell.prefix("source %s/env_PSV.cfg; " % SNAKEMAKE_DIR)

"""
Dependencies: should be taken care of by loading the env_PSV.cfg file
Requires a file called reads.orig.bam, ref.fasta should also be there, and it will take advantage 
of duplications.fasta, if it is there
"""

blasr = '~mchaisso/projects/AssemblyByPhasing/scripts/abp/bin/blasr'
blasrDir = '~mchaisso/projects/blasr-repo/blasr'
scriptsDir = '/net/eichler/vol5/home/mchaisso/projects/AssemblyByPhasing/scripts/abp'
#base2="/net/eichler/vol21/projects/bac_assembly/nobackups/scripts"
base2="/net/eichler/vol2/home/mvollger/projects/abp"

#settings
MINCOV = int(os.environ.get("MINCOV", "30"))
MAXCOV = int(os.environ.get("MAXCOV", "50"))
# only consider collapsed sites
MINTOTAL = int(os.environ.get("MINTOTAL", "100"))
#MINCOV=30
#MAXCOV=50
#MINTOTAL=100
print("MINCOV:{}\nMACCOV:{}\nMINTOTAL:{}".format(MINCOV, MAXCOV, MINTOTAL))

rule master:	
    input: 'PSV1_done'
    message: 'Running PSV1'

#
# This realigns reads with match/mismatch parameters that make it more
# likely to align a PSV as a mismatch rather than paired insertion and
# deletion.
#
rule preprocess_reads:
    input:
        'reads.orig.bam'
    output:
        'reads.bas.h5'
    shell:
        'samtools view -h {input} | {blasrDir}/pbihdfutils/bin/samtobas /dev/stdin {output}'

rule realign_reads:
    input:
        basreads='reads.bas.h5', 
        ref='ref.fasta'
    output:
        "reads.bam"
    shell: 
        '{blasr} {input.basreads} {input.ref} -sam  -mismatch 3 -insertion 9 -deletion 9 -nproc 4 -out /dev/stdout -minMapQV 30 -minAlignLength 2000 -preserveReadTitle | samtools view -bS - | samtools sort -T tmp -o {output}'

rule index_realigned_reads:
	input:
		"reads.bam"		
	output:
		"reads.bam.bai"	
	shell:
		'samtools index {input}'	


#
# Given the alignments, count the frequency of each base at every
# position, and estimate the SNVs from this frequency. 
#
#echo "{scriptsDir}/BamToSNVTable.sh reads.bam ref.fasta MINCOV" ??? what is this for ??? 
rule create_SNVtable_from_reads:
    input:
        reads="reads.bam",
        ref="ref.fasta"	
    output: 
        snv="assembly.consensus.fragments.snv",
        vcf="assembly.consensus.nucfreq.vcf",
    shell:
        '{scriptsDir}/BamToSNVTable.sh {input.reads} {input.ref} {MINCOV} {MINTOTAL}'


#
# Create a matrix with one row per read, and one column per PSV.  
#  . - read is ref base
#  1 - read is PSV
#  n - no signal (indel or read ended)
#
rule SNVtable_to_SNVmatrix:
    input:
        snv="assembly.consensus.fragments.snv"	
    output:
        mat="assembly.consensus.fragments.snv.mat",
        snvsPos="assembly.consensus.fragments.snv.pos"
    shell:
       '{scriptsDir}/FragmentSNVListToMatrix.py {input.snv} --named --pos {output.snvsPos} --mat {output.mat}'  


#
# create duplicaitons.fasta
#
GRCh38="/net/eichler/vol2/eee_shared/assemblies/GRCh38/GRCh38.no_alts.fasta"
sa="/net/eichler/vol2/eee_shared/assemblies/GRCh38/GRCh38.no_alts.fasta.sa"
rule duplicationsFasta1:
    input:
        dupref="ref.fasta"
    output:
        dupsam="ref.fasta.sam",
    threads: 8
    shell:
        '''
        blasr {input.dupref} {GRCh38} -nproc {threads} -sa {sa} -sam -bestn 30 -out {output.dupsam}
        '''

rule duplicationsFasta2:
    input:
        dupsam="ref.fasta.sam",
    output:
        duprgn="ref.fasta.rgn",
    shell:
        '''
        ~mchaisso/projects/mcutils/bin/samToBed {input.dupsam} --reportIdentity | awk '{{ if ($3-$2 >10000 && $9 > 0.85) print $1":"$2"-"$3 }}' > {output.duprgn}
        '''

rule duplicationsFasta3:
    input:
        duprgn="ref.fasta.rgn",
    output:
        dup="duplications.fasta",
    shell:
        '''
        samtools faidx {GRCh38} `cat {input.duprgn}` > {output.dup}
        '''



#
# Set up the ground truth if it exists.  Map the collapsed sequence
# back to the human genome, and output all th regions that are
# sufficiently similar to the collapse.
#
MAX_SEGDUPS="/net/eichler/vol5/home/mchaisso/projects/SegDups/GRCh38/grch38_superdups.max_identity.bed"
MEAN_SEGDUPS="/net/eichler/vol5/home/mchaisso/projects/SegDups/GRCh38/grch38_superdups.mean_identity.bed"

rule reads_to_fasta:
    input:
        "reads.bam"
    output:
        "reads.fasta"
    shell:
        '''samtools view {input} | awk '{{ print ">"$1; print $10; }}' > {output}'''

rule realignReads_to_Dups:
    input:
        "reads.fasta",
        "duplications.fasta"
    output:
        "reads.dups.m4"
    shell:
        '{blasr} {input} -m 4 -bestn 1 -preserveReadTitle -out {output} -nproc 4'

rule orderMatByalignments: #get more explanation
    input:
        "assembly.consensus.fragments.snv.mat",
        "reads.dups.m4"
    output:
        "assembly.consensus.fragments.snv.mat.categorized"
    shell:
        '{scriptsDir}/sorting/OrderMatByAlignments.py {input}  > {output}'

rule fastaToBed:
    input: "duplications.fasta"
    output: "duplications.bed"
    run:
        fasta = list(SeqIO.parse("duplications.fasta", "fasta"))
        rtn = ""
        for rec in fasta:
            temp = rec.id.split(":")
            chro = temp[0]
            nums = temp[1].split("-")
            start = nums[0]
            end   = nums[1]
            rtn += "{}\t{}\t{}\n".format(chro, start, end)
        f = open(output[0], "w+")
        f.write(rtn)
        f.close()

rule sortBed:
    input: "duplications.bed"
    output: "duplications.sorted.bed"
    shell:
        '''
        bedtools sort -i {input} > {output}
        '''

rule intersectBed:
    input: "duplications.sorted.bed"
    output:
        maxx = "dup_max.intersect",
        mean = "dup_mean.intersect"
    shell:
        '''
        bedtools intersect -a {input} -b {MAX_SEGDUPS}  -wa -wb > {output.maxx}
        bedtools intersect -a {input} -b {MEAN_SEGDUPS} -wa -wb > {output.mean}
        '''

#if os.path.exists("duplications.fasta"):
if True:
    x="this is a palcehold in canse I need to put back in the conditional duplicaitons.fasta"
else:
    rule ifNoDuplicationsFasta:
        input:
            "assembly.consensus.fragments.snv.mat"
        output:
            "assembly.consensus.fragments.snv.mat.categorized"
        shell:
            '''
            {base2}/categorize.sh {input} {output}
            '''
            #'''cat {input} | awk { print $1"\t"$2"\tall" } > {output}'''

    rule noIntersect:
        input: 
            "assembly.consensus.fragments.snv.mat.categorized"
        output:
            maxx = "dup_max.intersect",
            mean = "dup_mean.intersect"
        shell:
            '''
            > {output.maxx}
            > {output.mean}
            '''

#
# This finds PSVs that are connected by a sufficient number of
# sequences, and creates the PSV graph. This will have merged components.
#
rule createPSVgraph:
	input:
		matrix="assembly.consensus.fragments.snv.mat.categorized",
		vcf="assembly.consensus.nucfreq.vcf"
	output:
	        graph="mi.gml",
                adj="mi.adj"
	shell:
		'{scriptsDir}/PairedSNVs.py {input.matrix} --minCov {MINCOV} --maxCov {MAXCOV} --graph {output.graph} --adj {output.adj} --minNShared 5 --minLRT 1.5 --vcf {input.vcf}'			


#
# Run correlation clustering to try and spearate out the graph.  Swap
# is a 'simulated annealing' type parameters. The factor parameter
# controls for how many negative edges are added for every positive edge.
#
rule correlationClustering:
	input:
		graph="mi.gml",
	output:
		out="mi.cuts.gml",
		plt="mi.gml.png",
		sites="mi.gml.sites",
                cuts="mi.gml.cuts"
	shell:
		'{scriptsDir}/MinDisagreeCluster.py --graph {input.graph} --cuts {output.cuts} --sites {output.sites} --factor 2 --swap 1000 --plot {output.plt} --out {output.out}'

rule index_ref_fasta:
	input:
		"ref.fasta"
	output:
		"ref.fasta.fai"
	shell:
		'samtools faidx {input}'




#
# Correlation clustering defines a set of cuts that separate out
# connected components.  This takes the cuts and makes a vcf file for
# each component.
#
rule makeCutsInPSVgraph:
    input:
        refIdx="ref.fasta.fai",
        cuts="mi.gml.cuts",
        snvsPos="assembly.consensus.fragments.snv.pos",
        vcf="assembly.consensus.nucfreq.vcf"	
    output:
        vcfs=dynamic("group.{n}.vcf"), 
    shell:
        '{scriptsDir}/CutsToPhasedVCF.py {input.cuts} {input.snvsPos} {input.vcf} --minComponent 4 --summary mi.comps.txt --ref {input.refIdx} '


#
# dynamic wildcards are grabing this across multiple directoeis, and there is no way I can figure out to 
# stop it from doing that so I made a second file so I dont have to use dynamic, this second snakemake file
# must be called by a master script or seperatly.
#
rule startPSV2:
    input: 
        dynamic("group.{n}.vcf"),
        maxx = "dup_max.intersect",
        mean = "dup_mean.intersect"
    output: "PSV1_done"
    shell: "touch {output}"



