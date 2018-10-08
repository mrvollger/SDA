#!/usr/bin/env bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $DIR/env_python3.cfg

# snakemake paramenters
snakefile=$DIR/ProcessCollapsedAssembly.py
jobNum=300
waitTime=60

# QSUB parameters, these are only the defualts, they can be changed with params.sge_opts
logDir=logs
mkdir -p $logDir
E=$logDir'/snakejob_{rule}_{wildcards}_e'
O=$logDir'/snakejob_{rule}_{wildcards}_o'

# run snakemake
# the first set of line runs it on a sub grid cluster
if [ "sungrid" == "sungrid" ]; then 
	snakemake -p \
		-s $snakefile \
		--drmaa " -P eichlerlab \
			-q eichler-short.q \
			-l h_rt=24:00:00  \
			-V -cwd -e $E -o $O \
			-l mfree={params.mem} \
			-pe serial {params.cores} \
			-S /bin/bash" \
		--jobs $jobNum \
		--latency-wait $waitTime \
		$1 $2 $3  # just a way to pass aditional arguments to snakemake, like --unlock 

# params.{cores,mem} goes and looks into the snakemake for the mem and cluster params and passes them as additional arguments

else
	snakemake -p \
		-s $snakefile \
		--jobs $(nproc) \
		 $1 $2 $3
fi 



