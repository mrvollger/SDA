#!/bin/bash

out=$PWD/collisionOfReads.txt
> $out
for region in $(cat regions.txt); do
	cd $region
	overlapOfReadsCheck.py "group*/H2.WH.bam" >> $out
	#touch reads.orig.bam
	#cp ../coverage.json coverage.json
	#summary.py
	#bedForABP.py Mitchell_CHM1
	cd ..
done

