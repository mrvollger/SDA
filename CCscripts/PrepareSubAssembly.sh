#!/usr/bin/env bash
~/projects/blasr/cpp/pbihdfutils/bin/samtobas   group.2.sam reads.bas.h5
blasr reads.bas.h5 assembly.quiver.fasta -sam -minMapQV 30 -minAlignLength 1000  -out /dev/stdout -nproc 4 | samtools view -bS - | samtools sort - -T tmp -o reads.bam
samtools index reads.bam
ln -s ../ref.dups.fasta .
ln -s assembly.quiver.fasta ref.fasta
