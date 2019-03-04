import os

snake_dir = os.path.dirname(workflow.snakefile) + "/"
shell.executable("/bin/bash")
#shell.prefix("source %s/env_python2.cfg; set -eo pipefail; " % snake_dir)
shell.prefix("source %s/env_python2.cfg; " % snake_dir)

#
# script locations and configurations 
#
if(os.path.exists("abp.config.json")):
	configFile = "abp.config.json"
else:
	configFile = "sda.config.json"
configfile:
	configFile	

MINCOV=config["MINCOV"]
MAXCOV=config["MAXCOV"]
MINTOTAL=config["MINTOTAL"]

base = snake_dir + "scripts/"
scriptsDir = snake_dir + "CCscripts/"
python3 = snake_dir + "env_python3.cfg"
python2 = snake_dir + "env_python2.cfg"


print(snake_dir)
print("MINCOV:{}\nMAXCOV:{}\nMINTOTAL:{}".format(MINCOV, MAXCOV, MINTOTAL))



onsuccess:
	print("SDA1 FINISHED: CC finished and created group.*.vcf")
onerror:
	print("SDA1 FAILED: Unable to create groups of phased PSVs.")


rule all:	
	input:
		pdf='CC/mi.cuts.gml.pdf',
		done = "CC/vcfs.done",
		png="Coverage.png",
		png2="CovWithCC.png",

#
# This realigns reads with match/mismatch parameters that make it more
# likely to align a PSV as a mismatch rather than paired insertion and
# deletion.
#
ISPB=True
ISONT=False
maptype= " map-pb "
if("ont" in config):
	if(config["ont"].lower() in ["t", "true"] ):
		ISPB=False
		ISONT=True
		maptype= " map-ont "

MM2 = False
if("minimap" in config):
	if(config["minimap"].lower() in ["t", "true"] ):
		MM2=True

bandwidth = 50000
# set a minimum alignment lenght 
# the minmum alignment length is set to be larger than a full length LINE element ~6 kbp. 
minaln = 3000
if("minaln" in config):
	minaln = int(config["minaln"])
minscore = int(minaln) # minimap does not have a minaln setting so I approximatie by setting a minium dynamic progrtaming score

if(os.path.exists("reads.orig.bam") and ISPB and not MM2):

	rule realign_reads:
		input:
			reads = 'reads.orig.bam',
			ref='ref.fasta',
		output:
			temp("reads.bam")
		threads:
			8
		shell: """
source {python2}
echo "Running Blasr"
blasr {input.reads} {input.ref} \
	--clipping subread \
	--sam \
	--out - \
	--nproc {threads} --bestn 1 \
	--mismatch 3 --insertion 9 --deletion 9 \
	--minAlignLength {minaln} \
	--minMatch 13 | \
	 samtools view -bS -F 4 - | \
	 samtools sort -m 4G -T tmp -o {output}

samtools view {output} > tmp.sam

# remove output if it makes an empty bam
if [ -s tmp.sam ] ; then
	rm tmp.sam
else
	rm tmp.sam
	rm {output}*
fi

"""

elif(os.path.exists("reads.orig.bam") and (ISONT or MM2)):
	rule realign_reads_minimap:
		input:
			reads='reads.orig.bam', 
			ref='ref.fasta'
		output:
			temp("reads.bam")
		threads: 8
		shell:"""
source {python3}

samtools fastq {input.reads} | \
	minimap2 \
		-ax {maptype} \
		--eqx \
		-L \
		-t {threads} \
		-k 11 \
		-A 3 \
		-B 3 \
		-O 9 \
		-E 3 \
		-s {minscore} \
		-r {bandwidth} \
		-R '@RG\\tID:BLASR\\tSM:NO_CHIP_ID\\tPL:PACBIO' \
		ref.fasta /dev/stdin | \
		samtools view -bS -F 2308 - | \
		samtools sort -m 4G -T tmp -o {output}
# cs adds a tag the generates info about insertion and deletion events 
# removed flaggs 4 (unmapped) + 256 (secondary) + 2048 (chimeric)
# actually i think I will allow chimeric, (allows jumps of large gaps)
# -A matching score -B mismatching score -E gap extenshion penalty 

"""



else:
	print("NO INPUT READS!!!")
	exit()


rule index_realigned_reads:
	input:
		"reads.bam"
	output:
		bai = temp("reads.bam.bai"),
	shell:
		'samtools index {input}'	



#
# this adds an index to each read name, make it so it works with WhatsHap
# It also adss a tlen that is the length of the alignment to the reference
# this is to mimic what balsr would have done 
# in addition it adds a SM tag to each read, this is for WhatsHap
#
rule AddIndextoReads:
	input:
		reads = "reads.bam",
		bai = "reads.bam.bai",
	output:
		reads = "reads.sample.bam",
		bai = "reads.sample.bam.bai",
	shell:"""
source {python3}
{base}/changeBamName.py {input.reads} {output.reads}
samtools index {output.reads}
"""



#
#
# lets just look at teh het profile
#
rule hetProfile:
	input:
		reads="reads.sample.bam",
		bai = "reads.sample.bam.bai",
		ref="ref.fasta"	,
	output: 
		nucfreq="snvs/nofilter.consensus.nucfreq",
	shell:
		"""
		samtools mpileup -q 0 -Q 0 {input.reads} | \
				{base}/MpileupToFreq.py  /dev/stdin | \
				{base}/PrintHetFreq.py 0 \
				--maxCount 10000000000000000 \
				--minTotal 0 \
				> {output.nucfreq}

		"""

rule thresholdProfile:
    input:
        nucfreq="snvs/nofilter.consensus.nucfreq",
    output: 
        png="Coverage.png",
    shell:
        """
		source {python3}
		{base}autoThreshold.py --nucfreq {input.nucfreq} --plot {output.png} 
        """

#
# Given the alignments, count the frequency of each base at every
# position, and estimate the SNVs from this frequency. 
#
rule create_nucfreq_from_reads:
	input:
		reads="reads.sample.bam",
		bai = "reads.sample.bam.bai",
		ref="ref.fasta",
	output:
		nucfreq="snvs/assembly.consensus.nucfreq",
	shell:"""
echo "Sam to nucfreq"
samtools mpileup -q 0 -Q 0 {input.reads} | \
	{base}/MpileupToFreq.py  /dev/stdin | \
	{base}/PrintHetFreq.py {MINCOV} \
	--maxCount {MINCOV} \
	--minTotal {MINTOTAL} \
	> {output.nucfreq}


if [ -s {output.nucfreq} ]; then
	echo "Candidate PSVs found!"
else
	>&2 echo "ERROR: NO candidate PSVS were found. SDA cannot resolve this duplication."
	rm {output}
	exit 1
fi
"""


rule create_SNVtable_from_reads:
	input:
		reads="reads.sample.bam",
		ref="ref.fasta"	,
		nucfreq="snvs/assembly.consensus.nucfreq",
	output: 
		snv="snvs/assembly.consensus.fragments.snv",
		vcf="snvs/assembly.consensus.nucfreq.vcf",
		filt="snvs/assembly.consensus.nucfreq.filt",
		frag="snvs/assembly.consensus.fragments",
	shell:"""
echo "filter nucfreq"
# this actaully only really creates the frag file, all filtering should be done by PrintHetfreq
samtools view -h reads.sample.bam | \
		{base}/readToSNVList  \
		--nft {input.nucfreq} \
		--sam /dev/stdin \
		--ref {input.ref} \
		--minFraction 0.01 \
		--minCoverage {MINCOV} \
		--out {output.frag} \
		--nftOut {output.filt}

echo "filtered to vcf"
{scriptsDir}/FreqToSimpleVCF.py --freq {output.filt} \
		--ref {input.ref} \
		--out {output.vcf} 

echo "vcf to snv"
{scriptsDir}/FragmentsToSNVList.py  \
		--fragment \
		--frags {output.frag} \
		--vcf {output.vcf} \
		--out {output.snv}


"""


#
# Create a matrix with one row per read, and one column per PSV.  
#  . - read is ref base
#  1 - read is PSV
#  n - no signal (indel or read ended)
#
rule SNVtable_to_SNVmatrix:
    input:
        snv="snvs/assembly.consensus.fragments.snv"	
    output:
        mat="snvs/assembly.consensus.fragments.snv.mat",
        snvsPos="snvs/assembly.consensus.fragments.snv.pos"
    shell:"""
{scriptsDir}/FragmentSNVListToMatrix.py {input.snv} --named --pos {output.snvsPos} --mat {output.mat}
"""


#
# Set up the ground truth if it exists.  Map the collapsed sequence
#
if( os.path.exists("duplications.fasta") and os.path.getsize("duplications.fasta") > 0 ):
	rule depthOnDuplications:
		input:
			reads="reads.sample.bam",
			ref = ancient( "duplications.fasta" ),
		output:
			bam = "dup/reads.dup.bam",
			depth="dup/dup_depth.tsv",
		threads:8
		shell:
			"""
			source {python3}
				samtools fastq {input.reads} | \
					minimap2 \
						-ax {maptype} \
						--eqx \
						-k 11 \
						-t {threads} \
						{input.ref} /dev/stdin | \
						samtools view -bS -F 2308 - | \
						samtools sort -m 4G -T tmp -o {output}
			
			samtools index {output}
			samtools depth -aa {output.bam} > {output.depth}
			"""


rule addFakeCatagoryToMatrix:
	input:
		"snvs/assembly.consensus.fragments.snv.mat"
	output:
		"snvs/assembly.consensus.fragments.snv.mat.categorized"
	shell:
		"""
		# this adds a fake catigory on the end
		cat {input} | awk '{{ print $1"\t"$2"\tall"}}' > {output}
		"""

#
# This finds PSVs that are connected by a sufficient number of
# sequences, and creates the PSV graph. This will have merged components.
#
if("MINLRT" in config):
	MINLRT = float(config["MINLRT"])
else:
	MINLRT = 1.5
rule createPSVgraph:
	input:
		matrix="snvs/assembly.consensus.fragments.snv.mat.categorized",
        png="Coverage.png",
		vcf="snvs/assembly.consensus.nucfreq.vcf"
	output:
		graph="CC/mi.gml",
		adj="CC/mi.adj",
		mi="CC/mi.mi",
	shell:"""
{scriptsDir}/PairedSNVs.py {input.matrix} --minCov {MINCOV} --maxCov {MAXCOV} \
		--mi {output.mi} --graph {output.graph} --adj {output.adj} \
		--minNShared 5 --minLRT {MINLRT} --vcf {input.vcf}

"""


#
# generate repulsion edges 
#
rule GenerateRepulsion:
	input:
		graph="CC/mi.gml",
		mi="CC/mi.mi",
	output:
		rep = "CC/mi.repulsion",
	shell:"""
{base}/GenerateRepulsion.py --shared 5 --lrt {MINLRT} --max 3 --gml {input.graph} --mi {input.mi} --out {output.rep}
"""

runIters = True
runIters = False
posFile = "CC/mi.gml"
if(runIters and os.path.exists(posFile)):
	print("Going to export iterations of CC for debugging")
	#showIters = " --exportEachIteration --layout {} ".format(posFile)	
	showIters = " --exportEachIteration "
	convert = "convert extraCCplots/iteration.*.png extraCCplots/all_cc_iterations.pdf"
else:
	showIters = ""
	convert = ""
#
# Run correlation clustering to try and spearate out the graph.  Swap
# is a 'simulated annealing' type parameters. The factor parameter
# controls for how many negative edges are added for every positive edge.
#
rule correlationClustering:
	input:
		graph="CC/mi.gml",
		rep = "CC/mi.repulsion",
	output:
		out="CC/mi.cuts.gml",
		sites="CC/mi.gml.sites",
		cuts="CC/mi.gml.cuts",
	shell:"""	
mkdir -p extraCCplots
rm -rf extraCCplots/*

{scriptsDir}/MinDisagreeClusterByComponent.py  \
	--graph {input.graph} \
	--repulsion {input.rep} \
	--niter 10000 --swap 1000000 --factor 1 \
	--plotRepulsion {showIters} \
	--cuts {output.cuts} --sites {output.sites} --out {output.out}

if [ -s {output.sites} ]; then
	echo "CC groups created!"
else
	>&2 echo "ERROR: No CC groups were made. SDA found no clusters of PSVs and cannot resolve this duplication."
	rm {output}
	exit 1
fi


# it running all iteraitons convert them into a booklet
{convert}

"""



rule index_ref_fasta:
	input:
		"ref.fasta"
	output:
		"ref.fasta.fai"
	shell:
		'samtools faidx {input}'


#
# Correlation clustering defines a set of cuts that separate out
# connected components.  This takes the cuts and makes a vcf file for
# each component.
#
rule makeCutsInPSVgraph:
	input:
		refIdx="ref.fasta.fai",
		cuts="CC/mi.gml.cuts",
		snvsPos="snvs/assembly.consensus.fragments.snv.pos",
		vcf="snvs/assembly.consensus.nucfreq.vcf"	
	output:
		#vcfs=dynamic("group.{n}.vcf"), # this caused a snakemake effore often because dynamic is depricated 
		done = "CC/vcfs.done",
	params:
		comps="CC/mi.comps.txt", 
	shell:
		"""
		rm -rf group.*
		{scriptsDir}/CutsToPhasedVCF.py {input.cuts} {input.snvsPos} {input.vcf} \
				--minComponent 4 \
				--ref {input.refIdx} \
				--summary {params.comps}
		touch {output.done}
		"""


#
# makes a gephi version of the plot, 
#
rule gephi:
	input:
		cuts="CC/mi.cuts.gml"
	output:
		pdf="CC/mi.cuts.gml.pdf",
	run:
		shell("mkdir -p extraCCplots")
		shell("source {python3}; {base}/gephi/gephi.sh {input.cuts} mi.cuts.gml" )
		shell("mv mi.cuts.gml.pdf {output.pdf}")
		collapse = os.path.basename(os.getcwd()) + ".pdf"
		shell("cp " + output["pdf"] + " " + collapse)




#
# makes a modified coverage plot with CC groups drawn on top
#
rule DepthProfileWithCC:
	input:
		nucfreq="snvs/nofilter.consensus.nucfreq",
		sites="CC/mi.gml.sites",
	output: 
		png="CovWithCC.png",
	shell:"""
source {python3}
{base}autoThreshold.py --nucfreq {input.nucfreq} --psvsites {input.sites} --plot {output.png} 
"""




