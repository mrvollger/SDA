#!/bin/bash


MEC.py assembly.consensus.fragments.snv.mat mi.gml.cuts --sam reads.bam --out reads.cc.bam
numGroups=$(wc -l < mi.gml.cuts)
maxCount=$(expr $numGroups - 1 )

for i in $(seq 0 $maxCount); do
	set -x 
	samtools view -h reads.cc.bam | grep "^@\|cc:i:$i" > group.$i/H2.WH.sam 
	set +x
done 

samtools view -h reads.cc.bam | grep "^@\|cc:i:-1" > reads.not.partitioned.sam

