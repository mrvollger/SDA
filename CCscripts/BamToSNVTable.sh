#!/usr/bin/env bash 
set -v
if [ $# -lt 2 ]; then
		echo "Usage: BamToSNVTable.sh file.bam ref.fasta mincov "
		exit 0
fi
chrom=""

PBS=/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/
samtools mpileup -q0 -Q 0 $1 | $PBS/Phasing/MpileupToFreq.py  /dev/stdin | $PBS/Phasing/PrintHetFreq.py $3 --maxCount $3 --minTotal $4 > assembly.consensus"$chrom".nucfreq 


echo "samtools view -h $1 |  /net/eichler/vol5/home/mchaisso/projects/pbgreedyphase/readToSNVList  --nft assembly.consensus.nucfreq --sam /dev/stdin  --ref $2 --out assembly.consensus.fragments --minFraction 0.01 --minCoverage $3 --nftOut assembly.consensus"$chrom".nucfreq.filt"

samtools view -h $1 |  /net/eichler/vol5/home/mchaisso/projects/pbgreedyphase/readToSNVList  --nft assembly.consensus.nucfreq --sam /dev/stdin  --ref $2 --out assembly.consensus.fragments --minFraction 0.01 --minCoverage $3 --nftOut assembly.consensus"$chrom".nucfreq.filt


~mchaisso/projects/AssemblyByPhasing/scripts/abp/FreqToSimpleVCF.py --freq assembly.consensus"$chrom".nucfreq.filt --out assembly.consensus."$chrom"nucfreq.vcf --ref $2

~mchaisso/projects/AssemblyByPhasing/scripts/abp/FragmentsToSNVList.py --frags assembly.consensus"$chrom".fragments --vcf assembly.consensus.nucfreq.vcf --out assembly.consensus.fragments.snv --fragment
