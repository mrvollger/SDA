#!/bin/bash

function stats {
    echo "Number of Assemblies"
    ls 0*/g*/"$1"assembly.fasta | wc -l
    echo "Assemblies of size 0" 
    grep "^>" 0*/g*/"$1"assembly.fasta -c | awk -F ":" '{print $2}' | grep -c "[0]"
    echo "Assemblies with two or more contigs"
    grep "^>" 0*/g*/"$1"assembly.fasta -c | awk -F ":" '{print $2}' | grep -c "[2-9]"
    echo "Assemblies with one contig"
    grep "^>" 0*/g*/"$1"assembly.fasta -c | awk -F ":" '{print $2}' | grep -c "[1]"
}


function stats2 {
    echo "Number of Assemblies"
    ls */group*.vcf | wc -l
    ls */group*/H2.WH.sam | wc -l
    ls */group*/"$1"assembly.consensus.fasta | wc -l
    echo "Assemblies of size 0" 
    grep "^>" */group*/"$1"assembly.consensus.fasta -c | awk -F ":" '{print $2}' | grep -c "[0]"
    echo "Assemblies with two or more contigs"
    grep "^>" */group*/"$1"assembly.consensus.fasta -c | awk -F ":" '{print $2}' | grep -c "[2-9]"
    echo "Assemblies with one contig"
    grep "^>" */group*/"$1"assembly.consensus.fasta -c | awk -F ":" '{print $2}' | grep -c "[1]"
}

#echo "GROUP 1"
#stats Mark. 
#echo

echo "WHATSHAP"
stats2 WH.

#echo "GROUP 2"
#stats BAD.
#echo ""

echo "Running on: "
qstat -u mvollger | grep -c "Run"
qstat -u mvollger | tail -n 1
echo ""


