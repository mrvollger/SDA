#!/usr/bin/env bash 
set -v
if [ $# -lt 2 ]; then
		echo "Usage: BamToSNVTable.sh file.bam ref.fasta mincov "
		exit 0
fi


chrom=".$4"


samtools mpileup -q0 -Q 0 $1 | /net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/Phasing/MpileupToFreq.py  /dev/stdin | $PBS/Phasing/PrintHetFreq.py $3 --maxCount $3 > assembly.consensus"$chrom".nucfreq
