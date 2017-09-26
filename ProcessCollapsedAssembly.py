import os
import tempfile

#
# setup the env for each exacution 
#
SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)
shell.executable("/bin/bash")
#shell.prefix("source %s/env_PSV.cfg; set -eo pipefail; " % SNAKEMAKE_DIR)
shell.prefix("source %s/env_PSV.cfg; " % SNAKEMAKE_DIR)

#
# A little complicated to find the temp dir
#
SSD_TMP_DIR = "/data/scratch/ssd"
if "TMPDIR" in os.environ:
    TMPDIR = os.environ['TMPDIR']
elif "TMPDIR" in config:
    TMPDIR = config['TMPDIR']
elif os.path.exists(SSD_TMP_DIR):
    TMPDIR = SSD_TMP_DIR
else:
    TMPDIR = tempfile.gettempdir()

configfile: "abp_on_denovo.json"


#wildcard_constraints:
    #index="\d+"


ruleorder: BedToCoverage > MergeCoverage > all

rule all:
	input:
		asmSA=config["asm"] + ".sa",
		split=dynamic("reads.{index}.fofn"),
		align=dynamic("alignments/align.{index}.bam"),
		alignbed=dynamic("alignments/align.{index}.bam.bed"),
		regions="coverage/regions.bed",
		#coverage=dynamic("coverage/coverage.{index}.bed"),
		combined="coverage/coverage.bed",
	# not used, just a test
	'''
	params:
		sge_opts="-P eichlerlab -pe serial 1 -l mfree=1G -l h_rt=1:00:00",
	output:
		"done"
	shell:
		"""
		echo {input.coverage}
		touch {output}
		"""
	'''

#
# If a suffix arry does not yet exist for the assembly, build this.
#
rule MakeASM_SA:
    input:
        asm=config["asm"]
    output:
        asmsa=config["asm"] + ".sa"
    params:
        sge_opts="-P eichlerlab -pe serial 1 -l mfree=48G -l h_rt=6:00:00",
        blasr=config["blasr"]
    shell:
        """
        {params.blasr}/alignment/bin/sawriter {input.asm}
        """

rule SplitFOFN:
    input:
        fofn=config["reads"]
    output:
        fofnSplit=dynamic("reads.{index}.fofn")
    params:
        sge_opts="-P eichlerlab -pe serial 1 -l mfree=1G -l h_rt=1:00:00",
        baxPerJob=config["bax_per_job"]
    shell:
        """
        split --numeric --lines {params.baxPerJob} {input.fofn} reads.
        for f in `ls reads*`; do
        mv $f $f.fofn
        done
        """

#
#  For read depth, and other future steps, it is necessary to map reads back to the assembly.
#
rule MapReads:
    input:
        asm=config["asm"],
        asmsa=config["asm"] + ".sa",
        fofn="reads.{index}.fofn"
    output:
        align="alignments/align.{index}.bam"
    params:
        sge_opts="-P eichlerlab -pe serial 8 -l mfree=8G -l h_rt=24:00:00",
        blasr=config["blasr"]
    shell:
        """
        mkdir -p alignments
        {params.blasr}/alignment/bin/blasr {input.fofn} {input.asm} \
            -nproc 8 -out /dev/stdout -minMapQV 0 -bestn 2 -advanceExactMatches 15 -sam | \
            samtools view -bS - | \
            samtools sort -T $TMPDIR/blasr -@ 8 -m4G - -o {output.align}
        samtools index {output.align}
        """
        

rule BamToBed:
    input:
        bam="alignments/align.{index}.bam"
    output:
        bed="alignments/align.{index}.bam.bed"
    params:
        sge_opts="-P eichlerlab -pe serial 1 -l mfree=2G -l h_rt=24:00:00",
        samutils=config["mcutils"]
    shell:
        """
        samtools view {input.bam} | {params.samutils}/samToBed /dev/stdin --reportIdentity > {output.bed}
        """

#
# this rule creats a bed file that is incremented by 100 for every contig
# these will be the feautes upon which we calculate depth wtih bedtools
#
rule FaiToBed:
    input:
        asmfai=config["asmfai"],
    output:
        regions="coverage/regions.bed",
    params:
        sge_opts="-P eichlerlab -pe serial 1 -l mfree=2G -l h_rt=24:00:00",
    run:
        shell("mkdir -p coverage")
        fai = open(input["asmfai"])
        out = ""
        for line in fai:
            token = line.strip().split("\t")
            length = int(token[1])
            contig = token[0]
            cur = 0
            for nxt in range(100, length, 100):
                out += "{}\t{}\t{}\n".format(contig, cur, nxt-1)
                cur = nxt
            out += "{}\t{}\t{}\n".format(contig, cur, length-1)
            #print(out)
            #print("{}\t{}".format(token[0], token[1]))
        outfile = open(output["regions"], "w+")
        outfile.write(out)



rule BedToCoverage:
    input:
        regions="coverage/regions.bed",
        bed="alignments/align.{index}.bam.bed",
    output:
        coverage="coverage/coverage.{index}.bed",
    params:
        sge_opts="-P eichlerlab -pe serial 1 -l mfree=2G -l h_rt=24:00:00",
    shell:
        """
        bedtools coverage -a {input.regions} -b {input.bed} -mean > {output.coverage}
        #bedtools merge -i test.sort.bed -c 4 -o sum
        """

rule MergeCoverage:
	input:
		coverage=dynamic("coverage/coverage.{index}.bed"),
	output:
		combined="coverage/coverage.bed",
	params:
		sge_opts="-P eichlerlab -pe serial 1 -l mfree=2G -l h_rt=24:00:00",
	shell:
		"""
		echo {input.coverage}
		#bedtools merge -i test.sort.bed -c 4 -o sum
		"""



