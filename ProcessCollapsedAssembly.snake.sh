#!/usr/bin/env bash
module purge
. /etc/profile.d/modules.sh
module load modules modules-init modules-gs/prod modules-eichler
export PATH=/net/eichler/vol2/home/mvollger/projects/builds/anaconda/anaconda3/bin:$PATH
set -x

#
# snakemake paramenters
#
snakefile=/net/eichler/vol2/home/mvollger/projects/abp/ProcessCollapsedAssembly.py
jobNum=200
waitTime=60 # this really needs to be 60 on our cluster :(
retry=1 # numer of times to retry the pipeline if it failes
# I allow a retry becuase sometimes even the really long waittime is not enough,
# and the files are actaully there

#
# QSUB parameters, these are only the defualts, they can be changed with params.sge_opts
# Allow snakemake to make directories, I think it slows things down when I done with "waitTime"
#
logDir=logs
mkdir -p $logDir
E=$logDir'/snakejob_{rule}_{wildcards}_e'
O=$logDir'/snakejob_{rule}_{wildcards}_o'
ram=4G
defaultCores=1

#
# run snakemake
#
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
	--restart-times $retry 


#--dryrun \
