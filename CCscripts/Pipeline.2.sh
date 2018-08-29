#!/usr/bin/env bash

source ~/scripts/setup_falcon.sh
module load mpc/0.8.2; module load mpfr/3.1.0; module load gmp/5.0.2; module load gcc/latest
module load numpy/1.11.0
module load scipy


if [ ! -e reads.bam ]; then
		samtools view -h reads.orig.bam | ~/projects/blasr/cpp/pbihdfutils/bin/samtobas /dev/stdin reads.bas.h5
		blasr reads.bas.h5 ref.fasta -sam  -mismatch 3 -insertion 9 -deletion 9 -nproc 4 -out /dev/stdout -minMapQV 30 -minAlignLength 2000 -preserveReadTitle | samtools view -bS - | samtools sort -T tmp -o reads.bam
		samtools index reads.bam
fi


/net/eichler/vol5/home/mchaisso/projects/AssemblyByPhasing/scripts/abp/BamToSNVTable.sh reads.bam ref.fasta


/net/eichler/vol5/home/mchaisso/projects/AssemblyByPhasing/scripts/abp/KPartition.py assembly.consensus.fragments.snv --named --pos assembly.consensus.fragments.snv.pos > assembly.consensus.fragments.snv.mat  


if [ -f duplications.fasta ]; then
		samtools view reads.bam | awk '{ print ">"$1; print $10;}' > reads.fasta
		blasr reads.fasta duplications.fasta -m 4 -bestn 1 -preserveReadTitle -out reads.dups.m4 -nproc 4
 		/net/eichler/vol5/home/mchaisso/projects/AssemblyByPhasing/scripts/abp/sorting/OrderMatByAlignments.py assembly.consensus.fragments.snv.mat reads.dups.m4  > assembly.consensus.fragments.snv.mat.categorized
else 
		cat assembly.consensus.fragments.snv.mat | awk '{ print $1"\t"$2"\tall"}' > assembly.consensus.fragments.snv.mat.categorized
fi


/net/eichler/vol5/home/mchaisso/projects/AssemblyByPhasing/scripts/abp/PairedSNVs.py assembly.consensus.fragments.snv.mat.categorized  --minCov 15 --maxCov 80  --graph mi.gml --adj mi.adj --minNShared 5 --minLRT 3 --vcf assembly.consensus.nucfreq.vcf 


samtools faidx ref.fasta
refname=`awk '{ print $1;}' ref.fasta.fai`
reflen=`awk '{ print $2;}' ref.fasta.fai`

/net/eichler/vol5/home/mchaisso/projects/AssemblyByPhasing/scripts/abp/ComponentsToPhasedVCF.py mi.gml assembly.consensus.fragments.snv.pos assembly.consensus.nucfreq.vcf --minComponent 4 --summary mi.comps.txt  --ref $refname --reflen $reflen --subgraph subgraph

for vcf in `ls group*.vcf`; do
		g=${vcf%.*}
		mkdir -p $g
    /net/eichler/vol5/home/mchaisso/projects/AssemblyByPhasing/scripts/abp/PartitionReads.sh $vcf $g/group 2;
		pushd $g
		/net/eichler/vol5/home/mchaisso/projects/AssemblyByPhasing/scripts/abp/AssembleRegion.2.sh group.2.sam
		popd
done

# Combine everything NOT in the alt as ref.
grep "^#" group.0.vcf > alts.vcf
cat group.*.vcf | grep -v "^#" | sort -k2,2n >> alts.vcf

mkdir -p ref
/net/eichler/vol5/home/mchaisso/projects/AssemblyByPhasing/scripts/abp/PartitionReads.sh alts.vcf ref/group 1;
pushd ref
/net/eichler/vol5/home/mchaisso/projects/AssemblyByPhasing/scripts/abp/AssembleRegion.2.sh group.1.sam
popd


ls group.*/assembly.quiver.fasta ref/assembly.quiver.fasta > assemblies.fofn

~/projects/AssemblyByPhasing/scripts/utils/RenameAssemblies.py assemblies.fofn assemblies.fasta
