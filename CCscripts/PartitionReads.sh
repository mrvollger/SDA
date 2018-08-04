 #!#/usr/bin/env bash
#
# Partition reads from phased vcf.
#
if [ $# -lt 2 ]; then
		echo "Usage: PartitionReads.sh partition.vcf outdir haplotype"
		exit 1
fi

base=$(dirname $0)

echo ""
echo "Start of PartitionReads.sh"
echo "Arguments:"
echo $1 $2 $3 
echo ""

# this command needs a --rgn argument added
echo "samtools view -h reads.bam \
    | ~mchaisso/projects/pbgreedyphase/partitionByPhasedSNVs \
        --vcf $1 \
        --ref ref.fasta \
        --h1 $2.1.sam  --h2 $2.2.sam --sam /dev/stdin \
        --unassigned /dev/null \
				--phaseStats $2.stats \
				--block 3 \
				--minDifference 0"
echo ""

samtools view -h reads.bam \
    | ~mchaisso/projects/pbgreedyphase/partitionByPhasedSNVs \
        --vcf $1 \
        --ref ref.fasta \
        --h1 $2.1.sam  --h2 $2.2.sam --sam /dev/stdin \
        --unassigned /dev/null \
				--phaseStats $2.stats \
				--block 4 \
				--minGenotyped 2 \
				--minDifference 3

echo ""
echo grep -v \"^@\" $2.$3.sam  | awk '{ print ">"$1; print $10;}' > $2.$3.fasta
grep -v "^@" $2.$3.sam  | awk '{ print ">"$1; print $10;}' > $2.$3.fasta
echo ""
echo "End of PartitionReads.sh"
echo ""
