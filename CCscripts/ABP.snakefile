configfile: "abp.json"
import os
import subprocess

rule all:
    "assembly.consensus.nucfreq.vcf"


rule DuplicationVCF:
		input:
		    "reads.bam"
		output:
		    "assembly.consensus.nucfreq.vcf"
		shell:
		    "/net/eichler/vol5/home/mchaisso/projects/AssemblyByPhasing/scripts/abp/BamToSNVTable.sh reads.bam ref.fasta"

rule SNVMat:
    input:
		    "assembly.consensus.nucfreq.vcf"
		output:
		    "assembly.consensus.fragments.snv.mat"
    shell:
        "/net/eichler/vol5/home/mchaisso/projects/AssemblyByPhasing/scripts/abp/KPartition.py assembly.consensus.fragments.snv --named --pos assembly.consensus.fragments.snv.pos --out assembly.consensus.fragments.snv.mat"

rule FASTAReads:
    input:
		    "reads.bam"
		output:
		    "reads.fasta"
		shell:
        samtools view reads.bam | awk '{ print ">"$1; print $10;}' > reads.fasta

rule CategorizedReads:
    input:
		    "duplications.fasta", "reads.fasta"
		output:
		    "reads.dups.m4"
		shell:
		    "blasr reads.fasta duplications.fasta -m 4 -bestn 1 -preserveReadTitle -out reads.dups.m4 -nproc 4"

rule CategorizedMat:
    input:
        "assembly.consensus.fragments.snv.mat", "reads.dups.m4"


    
