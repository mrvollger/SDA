#!/bin/bash
#$ -P eichlerlab
#$ -S /bin/bash -V
#$ -l mfree=10G
#$ -l h_rt=06:00:00
#$ -pe serial 1
#$ -cwd
#$ -q eichler-short.q


module load latex2rtf/2.3.2a


pdflatex --draftmode $1
pdflatex $1

