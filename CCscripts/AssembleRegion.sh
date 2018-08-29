#!/usr/bin/env bash
grep -v "^@" $1 | awk '{ print ">"$1;print $10;}' > reads.fasta
# There seems to be an error associated with the header of reads.fasta
# Start Changes by Mitchell on 05/08/2017
rm -f reads.fasta
#sed -i -e 's/_//g' reads.fasta
#sed -i -e 's/\///g' reads.fasta
# End changes by Mitchell on 05/08/2017

make -f /net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/local_assembly/CanuSamAssembly.mak assembly.consensus.fasta SAM=$1

if [ ! -s assembly.fasta ]; then
		rm -f assembly.fasta
		rm -f reads.fasta
		make -f /net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/local_assembly/FalconSamAssembly.mak assembly.consensus.fasta SAM=$1

fi

if [ ! -s assembly.fasta ]; then
		#
		# If the assembly didn't work, just try and call a quiver consensus.  This may not correctly
		# reflect structural variation, but it can add a few regions as consensus.
		#
		QuiverDir=/net/eichler/vol5/home/mchaisso/software/quiver
		unset QRSH_COMMAND && source $QuiverDir/setup_quiver.sh
		
		/net/eichler/vol5/home/mchaisso/projects/blasr/cpp/pbihdfutils/bin/samtobas $1 reads.bas.h5 -defaultToP6
		blasr reads.bas.h5 ../ref.fasta -sam -bestn 1 -nproc 4 -out /dev/stdout -clipping soft | samtools view -bS - | samtools sort -T tmp -o reads.to_ref.bam
		samtools index reads.to_ref.bam
		covStart=`samtools depth reads.to_ref.bam | head -1 | cut -f 1`
		covEnd=`samtools depth reads.to_ref.bam | tail -1 | cut -f 1`
		refName=`head -1 ../ref.fasta | tr -d ">\n"`
		/net/eichler/vol5/home/mchaisso/software/quiver/bin/quiver  -j8 --minCoverage 1 --noEvidenceConsensusCall nocall --referenceFilename ../ref.fasta reads.to_ref.bam -o assembly.quiver.orig.fasta
		/net/eichler/vol5/home/mchaisso/projects/AssemblyByPhasing/scripts/abp/RemoveFlankingNs.py assembly.quiver.orig.fasta assembly.quiver-only.fasta

fi
		
		



