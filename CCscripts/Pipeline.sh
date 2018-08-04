#!/usr/bin/env bash

source ~mchaisso/scripts/setup_abp
#source config.sh

MINCOV=30
MAXCOV=50
# only consider collapsed sites
MINTOTAL=100

while [[ $# -gt 1 ]]; do
		key="$1"
		case $key in
				-n|--minCov)
						MINCOV="$2"
						shift
						;;
			  -m|--maxCov)
						MAXCOV="$2"
						shift
						;;
				-h)
						echo "Usage: CreateBrowserShots.sh gaps.bed aln_dir dest_dir "
						echo "  -n, --minCov val Minimal read count for an alternate allele to be considered a PSV."
						echo "  -x, --maxCov val Maximal read count for an alternate allele to be considered a PSV."
						exit 1						
						;;
				*)
						;;
		esac
		shift
		
done
base=$(dirname $0)

#
# This realigns reads with match/mismatch parameters that make it more
# likely to align a PSV as a mismatch rather than paired insertion and
# deletion.
#
if [ ! -e reads.bam ]; then
		samtools view -h reads.orig.bam | $base/blasr/pbihdfutils/bin/samtobas /dev/stdin reads.bas.h5
		$base/bin/blasr reads.bas.h5 ref.fasta -sam  -mismatch 3 -insertion 9 -deletion 9 -nproc 4 -out /dev/stdout -minMapQV 30 -minAlignLength 2000 -preserveReadTitle | samtools view -bS - | samtools sort -T tmp -o reads.bam
		samtools index reads.bam
fi

#
# Given the alignments, count the frequency of each base at every
# position, and estimate the SNVs from this frequency. 
#
echo "$base/BamToSNVTable.sh reads.bam ref.fasta $MINCOV"
$base/BamToSNVTable.sh reads.bam ref.fasta $MINCOV $MINTOTAL

#
# Create a matrix with one row per read, and one column per PSV.  
#  . - read is ref base
#  1 - read is PSV
#  n - no signal (indel or read ended)
#

$base/FragmentSNVListToMatrix.py assembly.consensus.fragments.snv --named --pos assembly.consensus.fragments.snv.pos --mat assembly.consensus.fragments.snv.mat  

#
# Set up the ground truth if it exists.  Map the collapsed sequence
# back to the human genome, and output all th regions that are
# sufficiently similar to the collapse.
#
if [ -f duplications.fasta ]; then
		samtools view reads.bam | awk '{ print ">"$1; print $10;}' > reads.fasta
		# CHANGE 07/24/2017 Mitchell Vollger, I am 95 % sure this dir is wrong. changed $base/alignment/bin/blasr to $base/blasr/alignment/bin/blasr
        # Change continued, this version of blasr jsut breaks for me, so I am just using the blasr command that is at line 43.        # the oroginal line is below commented out  
        #$base/alignment/bin/blasr reads.fasta duplications.fasta -m 4 -bestn 1 -preserveReadTitle -out reads.dups.m4 -nproc 4
        $base/bin/blasr reads.fasta duplications.fasta -m 4 -bestn 1 -preserveReadTitle -out reads.dups.m4 -nproc 4
        # END CHANGE 07/24/2017 Mitchell Vollger
 		$base/sorting/OrderMatByAlignments.py assembly.consensus.fragments.snv.mat reads.dups.m4  > assembly.consensus.fragments.snv.mat.categorized
else 
		cat assembly.consensus.fragments.snv.mat | awk '{ print $1"\t"$2"\tall"}' > assembly.consensus.fragments.snv.mat.categorized
fi

#
# This finds PSVs that are connected by a sufficient number of
# sequences, and creates the PSV graph.  This will have merged components.
#
$base/PairedSNVs.py assembly.consensus.fragments.snv.mat.categorized  --minCov $MINCOV --maxCov $MAXCOV  --graph mi.gml --adj mi.adj --minNShared 5 --minLRT 1.5 --vcf assembly.consensus.nucfreq.vcf

#
# Run correlation clustering to try and spearate out the graph.  Swap
# is a 'simulated annealing' type parameters. The factor parameter
# controls for how many negative edges are added for every positive edge.
#

$base/MinDisagreeClusterByComponent.py --graph mi.gml --cuts mi.gml.cuts --sites mi.gml.sites --factor 2 --out mi.cuts.gml --swap 1000 --plot mi.gml.png --out mi.cuts.gml --plotRepulsion

samtools faidx ref.fasta

#
# Correlation clustering defines a set of cuts that separate out
# connected components.  This takes the cuts and makes a vcf file for
# each component.
#

$base/CutsToPhasedVCF.py mi.gml.cuts assembly.consensus.fragments.snv.pos assembly.consensus.nucfreq.vcf --minComponent 4 --summary mi.comps.txt  --ref ref.fasta.fai


#
# Now partition reads using each vcf, and assemble the reads in each
# partition.
#
for vcf in `ls group*.vcf`; do
		g=${vcf%.*}
		mkdir -p $g
		n=`echo $vcf | cut -d. -f 2`
		$base/PartitionReads.sh $vcf $g/group 2;
		pushd $g
		# CHANGE 07/24/2017, Mitchell Vollger changed group.2.sam to group.1.sam
        $base/AssembleRegion.sh group.1.sam
        # END CHANGE 07/24/2017 Mitchell Vollger
		popd
done

# CHANGE 07/24/2017 Mitchell Vollger
# the first one is no longer called assembly.quiver it is assembly.consensous 
ls group.*/assembly.consensus.fasta ref/assembly.quiver.fasta > assemblies.fofn
# END CAHNGE

# CHANGE -7/24/2017 Mitchell Vollger
# this script requires that --fofn and --out are included and they were not originally 
~mchaisso/projects/AssemblyByPhasing/scripts/utils/RenameAssemblies.py --fofn assemblies.fofn --out assemblies.fasta
# END CHANGE


