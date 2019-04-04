#!/usr/bin/env bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $DIR/env_python3.cfg

# snakemake paramenters
snakefile=$DIR/Bax2Bam.snakemake.py
jobNum=100
waitTime=60 # this really needs to be 60 on our cluster :(

# QSUB parameters, these are only the defualts, they can be changed with params.sge_opts
logDir=logs
mkdir -p $logDir
E=$logDir'/snakejob_{rule}_{wildcards}_e'
O=$logDir'/snakejob_{rule}_{wildcards}_o'
ram=8G
defaultCores=1

#
# run snakemake
# the first set of line runs it on a sun grid cluster
#
if [ "sungrid" == "sungrid" ]; then 
	snakemake -p \
		-s $snakefile \
		--drmaa " -P eichlerlab \
			-q eichler-short.q \
			-l h_rt=24:00:00  \
			-l mfree=$ram \
			-V -cwd -e $E -o $O \
			{params.cluster} \
			-S /bin/bash" \
		--jobs $jobNum \
		--latency-wait $waitTime \
		$1 $2 $3  # just a way to pass aditional arguments to snakemake, like --unlock 
else
	snakemake -p \
		-s $snakefile \
		--jobs $(nproc) \
		 $1 $2 $3
fi 



