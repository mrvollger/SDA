SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)
shell.executable("/bin/bash")
shell.prefix("source %s/env_PSV.cfg; " % SNAKEMAKE_DIR)

#parser = argparse.ArgumentParser(description="")
#parser.add_argument("--cuts", help="this shoudl be mi.gml.cuts", default="mi.gml.cuts" )
#parser.add_argument("--group", help="", type=int, default=None )
#parser.add_argument("--mat", help="this shoudl be the mat", default="assembly.consensus.fragments.snv.mat" )
#parser.add_argument("--sam", help="this shoudl be the sam", default="reads.bam" )
#parser.add_argument("--fasta", help="this shoudl be the sam", default="reads.fasta" )
#args = parser.parse_args()
import pysam 
from Bio import SeqIO
import pandas as pd

cutsfile = open( "mi.gml.cuts" ).readlines()
cuts = {}
for idx, cut in enumerate( cutsfile ):
	cutl =  cut.split()
	cut = [ int(x) for x in cutl ]
	cuts[idx] = cut

groups = cuts.keys()

def readsThatAgree( mat, cut ):
	mat = open(mat).readlines()
	reads = {}
	for read in mat:
		readpos = [ read[x] for x in cut ]
		# make sure there is at least one PSV, and that most of the PSVs are of the right group
		if('1' in readpos and ( readpos.count("1") > readpos.count(".") ) ):
			vector, key = read.split()
			reads[key] = readpos
	return(reads)


def readSam(reads, group, sam ):
	samfile = pysam.AlignmentFile( sam )
	samreads = {}
	for read in samfile.fetch(until_eof=True):
		key = read.query_name
		if(key in reads and read.infer_query_length() > 15000 ):
			samreads[key] = read
	
	for key in samreads:
		#print(key, samreads[key])
		print(samreads[key].infer_query_length())


def readFasta(reads, fasta, outname):
	fasta = list( SeqIO.parse(fasta, "fasta" ) )
	outreads = []	
	for read in fasta:
		key = read.id
		if(key in reads and len(read.seq) > 5000 ):
			outreads.append(read)
	
	SeqIO.write(outreads, outname, "fasta" ) 


#
# start snake
#
rule all:
	input:
		"done.txt"


rule getReads:
	input:
		cuts = "mi.gml.cuts",
		mat = "assembly.consensus.fragments.snv.mat",
		reads = "reads.fasta",
	output:
		reads = expand("extendEdges/{id}.fasta", id=groups ),
	run:
		for group in cuts:
			reads = readsThatAgree(input["mat"], cuts[group])
			outFile = "extendEdges/{}.fasta".format(group)
			readFasta(reads, input["reads"], outFile)



rule mapReads:
	input:
		reads = "extendEdges/{id}.fasta",
		asm = "group.{id}/WH.assembly.consensus.fasta",
	output:
		blasr = "extendEdges/{id}.blasr"
	shell:
		"""
		blasr -bestn 1 -m 5 {input.reads} {input.asm} > {output.blasr}
		"""

rule overhangReads:
	input:
		reads = "extendEdges/{id}.fasta",
		blasr = "extendEdges/{id}.blasr"
	output:
		reads = "extendEdges/{id}.overhang.fasta", 
	run:
		cnames = "qName qLength qStart qEnd qStrand tName tLength tStart tEnd tStrand score numMatch numMismatch numIns numDel mapQV qAlignedSeq matchPattern tAlignedSeq".split()

		df = pd.read_csv(input["blasr"], sep=" ", header = None)
		df.columns = cnames 
		toKeep = {}
		for idx, read in df.iterrows():
			thresh = 2000
			tTail = read["tLength"] - read["tEnd"]
			qname = "/".join( read["qName"].split("/")[:-1] )
			
			if(read["tStart"] < thresh ): # check if it overlaps the front
				head = read["qStart"] - read["tStart"]	# amount of sequence that extends past contig
				if(head > thresh): # the overlpa past the end of the contig is large enough
					toKeep[qname] = (0, head)
					print(read)
					continue 	
			
			if( tTail < thresh ):
				tail = read["qLength"] - read["qEnd"] - tTail
				if(tail > thresh):
					toKeep[qname] = ( read["qEnd"] + tTail, read["qLength"]   )
					print(read)


			if(  (read["tStart"] < thresh) or (tTail < thresh ) and False):
				if( read["qStart"] > read["tStart"] ):
					toKeep[qname] = (0, read["qStart"] )
					#print(read["qStart"], read["tStart"])
				elif( (read["qLength"]- read["qEnd"]) > tTail  ):
					toKeep[qname]= (read["qEnd"], read["qLength"] )
					#print( read["qLength"] - read["qEnd"], tTail  )
		print(toKeep)
		reads = list(SeqIO.parse(input["reads"], "fasta"))
		outReads = []
		for read in reads:
			if(read.name in toKeep.keys()):
				start, end = toKeep[read.name]
				read.seq = read.seq[start:end]
				outReads.append(read)
		SeqIO.write(outReads, output["reads"], "fasta" )


rule Next:
	input:
		expand = expand("extendEdges/{id}.overhang.fasta", id=groups ),
	output:
		done = "done.txt"
	shell:
		"""
		# touch {output.done}
		"""



