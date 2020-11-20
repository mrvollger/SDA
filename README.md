# Segmental Duplication Assembler (SDA)


# Download: #
```
git clone --recurse-submodules https://github.com/mrvollger/SDA.git
```

# Install: #
The requirements for SDA are taken care of by two custom conda environments (sda-python-2 and sda-python-3).
In order to run SDA, you must already have anaconda 3 installed on your system and you must be able to create conda environments. 

Once that is done, create `env_sda.sh` so that it adds conda to your path.  
Additionally, RepeatMasker, gcc, and cmake are not installed for you via conda and must be loaded in `env_sda.sh`.
Finally, if for example RepeatMasker relies on a perl module you must also add that to `env_sda.sh`
	Notes:
		1) To install `Racon` cmake must be avalible.
		2) To run `SDA denovo` there must be a version of `RepeatMasker` with the `ncbi` engine avalible.
		3) We suggest using `gcc=6.4.0`.
 
Here is an example of what your `env_sda.sh` might look like:
```
#!/bin/bash

unset PYTHONPATH
module purge
. /etc/profile.d/modules.sh
module load modules modules-init modules-gs/prod modules-eichler
module load gcc/6.4.0
module load cmake/3.14.3
module load perl/5.14.2
module load RepeatMasker/4.1.0

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

source /net/eichler/vol26/projects/sda_assemblies/nobackups/software/miniconda3/etc/profile.d/conda.sh
conda activate
```

Once `env_sda.sh` has been created the `Makefile` can be run with:
```
make
```
Several people have run into an error with readToSNVList if it is complied with gcc 8.x so please complie with gcc 6.4.0. 


# Run: #

```
./SDA  -h
usage: SDA <command> [<args>]

Segmental Duplication Assembler (SDA) commands options:
        denovo          Run SDA on a denovo assembly.
        collapse        Run SDA on a specific collapsed region.

Please cite:
        Vollger MR, et al. Long-read sequence and assembly of segmental duplications.
        Nat Methods. 2019 Jan;16(1):88-94. doi: 10.1038/s41592-018-0236-3.

        PMCID: PMC6382464.

positional arguments:
  {denovo,collapse}  denovo: run SDA on a denovo assembly while distributing
                     to a cluster. collapse: run SDA on a single collapsed
                     region.

optional arguments:
  -h, --help         show this help message and exit
```

## Run "SDA collapse" ##
```
./SDA collapse -h
usage: SDA [-h] [--coverage COVERAGE] [--reads READS] [--ref REF]
           [--platform {subread,ccs,ont,SUBREAD,CCS,ONT}] [-t THREADS]
           [-d DIR] [-p PREFIX] [--minaln MINALN] [--bandwidth BANDWIDTH]
           [--iterations ITERATIONS] [--assemblers ASSEMBLERS] [--lrt LRT]
           [--minNumShared MINNUMSHARED] [--maxPosRep MAXPOSREP]
           [--minCutSize MINCUTSIZE] [--minCutLen MINCUTLEN] [--debug]

SDA collapse

optional arguments:
  -h, --help            show this help message and exit
  --coverage COVERAGE   The average aligned read depth of the genome (default:
                        None)
  --reads READS         file with reads in it (default: reads.orig.bam)
  --ref REF             reference fasta file (default: ref.fasta)
  --platform {subread,ccs,ont,SUBREAD,CCS,ONT}
                        type of long read. (default: subread)
  -t THREADS, --threads THREADS
                        Threads for the snakemake (collapse), jobs for
                        snakemake (denovo) (default: 8)
  -d DIR, --dir DIR     directory for output files (default: sda_out)
  -p PREFIX, --prefix PREFIX
                        prefix for output files (default: sda)
  --minaln MINALN       Minimum alignment length (default: 3000)
  --bandwidth BANDWIDTH
                        bandwidth used in alignment (default: 50000)
  --iterations ITERATIONS
                        Number of times to run CC (default: 10)
  --assemblers ASSEMBLERS
                        Which assemblers to use for local assembly. canu or
                        wtdbg2 or both with commas and no spaces. Assemblies
                        are reported from only the first assembler unless it
                        fails, in which case the second assembler is used and
                        so on. (default: canu,wtdbg2)
  --lrt LRT             Required log likelihood ratio of reads with PSV vs
                        reads without (default: 1.5)
  --minNumShared MINNUMSHARED
                        The minimum number of reads that must span between two
                        PSVs (default: 5)
  --maxPosRep MAXPOSREP
                        The maximum number of reads that can link two PSVs and
                        still mark an edge as negative (default: 3)
  --minCutSize MINCUTSIZE
                        The minimum number of PSVs in a cluster (default: 4)
  --minCutLen MINCUTLEN
                        The minimum distance spanned by the PSVs in a cluster
                        (default: 9000)
  --debug               If set the temporary files are not deleted. (default:
                        False)
```

## Run "SDA denovo" ##
```
./SDA denovo -h
usage: 
		SDA denovo --input <input.(fofn|bam)> [<args>]
		For allowing cluster submission please add --cluster or --drmaa, these areguments are passed directly to snakemake. 
		The cluster/drmaa string must include these string: {threads} and {resources.mem}G . Below is an example using drmaa and SGE:
	
	--drmaa " -l mfree={resources.mem}G -pe serial {threads} -l h_rt=128:00:00 -V -cwd -S /bin/bash " 

SDA denovo

optional arguments:
  -h, --help            show this help message and exit
  --input INPUT         file of file names to align to the genome, or a bam
                        containg reads aligned to the genome (softclippled
                        recomended). (default: None)
  --species SPECIES     species name of data base for repeat makser (default:
                        human)
  --cluster CLUSTER     cluster configuration line for snakemake (default:
                        None)
  --drmaa DRMAA         drmaa configuration line for snakemake (default: None)
  --ref REF             reference fasta file (default: ref.fasta)
  --platform {subread,ccs,ont,SUBREAD,CCS,ONT}
                        type of long read. (default: subread)
  -t THREADS, --threads THREADS
                        Threads for the snakemake (collapse), jobs for
                        snakemake (denovo) (default: 8)
  -d DIR, --dir DIR     directory for output files (default: sda_out)
  -p PREFIX, --prefix PREFIX
                        prefix for output files (default: sda)
  --minaln MINALN       Minimum alignment length (default: 3000)
  --bandwidth BANDWIDTH
                        bandwidth used in alignment (default: 50000)
  --iterations ITERATIONS
                        Number of times to run CC (default: 10)
  --assemblers ASSEMBLERS
                        Which assemblers to use for local assembly. canu or
                        wtdbg2 or both with commas and no spaces. Assemblies
                        are reported from only the first assembler unless it
                        fails, in which case the second assembler is used and
                        so on. (default: canu,wtdbg2)
  --lrt LRT             Required log likelihood ratio of reads with PSV vs
                        reads without (default: 1.5)
  --minNumShared MINNUMSHARED
                        The minimum number of reads that must span between two
                        PSVs (default: 5)
  --maxPosRep MAXPOSREP
                        The maximum number of reads that can link two PSVs and
                        still mark an edge as negative (default: 3)
  --minCutSize MINCUTSIZE
                        The minimum number of PSVs in a cluster (default: 4)
  --minCutLen MINCUTLEN
                        The minimum distance spanned by the PSVs in a cluster
                        (default: 9000)
  --debug               If set the temporary files are not deleted. (default:
                        False)

```

There is a test case for this snakemake but it is rather large, 20Gb of data. To download it execute the following:
```
make TestCases/GenomeTest/ref.fasta
```


## Run SDA to identify regions of collapse ##
If you only need to run SDA for collapse analysis you can do a simplified install by first creating `env_sda.sh` which must add RepeatMasker (version 4.1.0) to your path and activate conda if it is not already in your default env, e.g.: 
```
PATH=$PATH:/path/to/your/repeat/masker/bin
```

You can then use make to install the other dependencies:
```
make collapse 
```

Run SDA only to the point of identifying regions of collapse (example). 
```
./SDA denovo \
        --drmaa " -l mfree={resources.mem}G -pe serial {threads} -l h_rt=128:00:00 -V -cwd -S /bin/bash " \
        --threads 400 \
        --platform ccs \
        --fofn /net/eichler/vol27/projects/sequence_data/nobackups/human/CHM13/PacBioHiFi/20_kbp_insert_hifi_beta/fastq.fofn \
        --ref ref.fasta \
		--pre sda  \
        sda_out/coverage/sda.collapses.bed
```



There are two collapse files output by SDA both with the following format:
```
contig	start	end	mean_coverage	median_coverage	#_of_bases_with_common_repeat_elements	length
```
The file named like `*.collapses.with.cm.bed` has all collapses regardless of common repeat elements. The file named `*.collapses.bed` will only have collapses with less than 75% common repeat element. 


The following commands makes a table with collapsed and expanded bases counted. 
```
./scripts/count_collapse.py --coverage 25 sda_out/coverage/*collapses.*.bed
```
### Note on HiFi reads ###
While SDA does now have presets for assembling with HiFi data it has not been extensively optimized for HiFi. I recommend using HiCanu (doi:10.1101/2020.03.14.992248) for assembling highly identical duplicates with high accuracy long reads. 

### Please Cite ###
Vollger, M. R. et al. Long-read sequence and assembly of segmental duplications. Nature Methods (2018). doi:10.1038/s41592-018-0236-3


