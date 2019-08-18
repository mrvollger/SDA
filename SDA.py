#!/usr/bin/env python
import os
import sys
import argparse
import subprocess 

# DEFAULTS 
minNumShared = 5
maxPosRep = 3
lrt = 1.5
minCutSize = 4
minCutLen = 9000

parser = argparse.ArgumentParser(description="""Segmental Duplication Assembler (SDA). \t 
Please cite: \t Vollger MR, et al. Long-read sequence and assembly of segmental duplications. Nat Methods. 2019 Jan;16(1):88-94. doi: 10.1038/s41592-018-0236-3. PMCID: PMC6382464.
""",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--coverage", help="The average aligned read depth of the genome", type=int)
parser.add_argument("--platform", help="type of long read.", default="subread", choices=["subread", "ccs", "ont"])
parser.add_argument("-t", "--threads", help="Threads for the snakemake", default=8, type=int)
parser.add_argument("--reads", help="file with reads in it", default="reads.orig.bam" )
parser.add_argument("--ref", help="reference fasta file", default="ref.fasta")
parser.add_argument("-d", "--dir", help="directory for output files", default="sda_out")
parser.add_argument("-p","--prefix", help="prefix for output files", default="sda")
parser.add_argument("--minaln", help="Minimum alignment length", default=3000 , type=int)
parser.add_argument("--bandwidth", help="bandwidth used in alignment", default=50000 , type=int)
parser.add_argument("--iterations", help="Number of times to run CC", default=10 , type=int)
parser.add_argument("--assemblers", help="""
	Which assemblers to use for local assembly. canu or wtdbg2 or both with commas and no spaces. 
	Assemblies are reported from only the first assembler unless it fails,
	in which case the second assembler is used and so on.""" ,
	default="canu,wtdbg2")
parser.add_argument("--lrt", help="Required log likelihood ratio of reads with PSV vs reads without", default=lrt , type=float)
parser.add_argument("--minNumShared", help="The minimum number of reads that must span between two PSVs", default=minNumShared , type=int)
parser.add_argument("--maxPosRep", help="The maximum number of reads that can link two PSVs and still mark an edge as negative", 
		default=maxPosRep, type=int)
parser.add_argument("--minCutSize", help="The minimum number of PSVs in a cluster", default=minCutSize, type=int)
parser.add_argument("--minCutLen", help="The minimum distance spanned by the PSVs in a cluster", default=minCutLen, type=int)
args, unknown = parser.parse_known_args()

if( not args.coverage ):
	parser.error('No coverage specified, add --coverage {\d+}')
if(not os.path.exists(args.ref + ".fai")):
	parser.error(f'{args.ref} must be indexed, try samtools faidx {args.ref}')


if(args.platform == "ccs"):
	#CCS DEFAULTS, if unchanged by used 
	if(args.minNumShared == minNumShared):
		minNumShared = 2
	if(args.maxPosRep == maxPosRep ):
		maxPosRep = 1
	if(args.lrt == lrt):
		lrt = 1.0
	if(args.minCutSize == minCutSize):
		minCutSize = 2
	if(args.minCutLen == minCutLen):
		minCutLen = 100 


# unknown paramters will be passed to snakemake 
snakeargs = " ".join(unknown) 
if(len(unknown) > 0):
	sys.stderr.write(f"Extra arguments passed to snakemake: {snakeargs}\n\n")


DIR = os.path.dirname( os.path.realpath(__file__) )
CWD = os.getcwd()

cmd=f'''cd {DIR} && source {DIR}/env_python3.sh && cd {CWD} \
&& snakemake \
-p -j {args.threads} -s {DIR}/SDA.smk {snakeargs} \
--config \
'''

for arg in vars(args):
	if(arg not in ["threads"] ):
		val = getattr(args, arg)
		cmd += f'{arg}={val} '

#sys.stderr.write(f"Snakemake command:\n{cmd}\n\n")
subprocess.call(cmd, shell=True)

