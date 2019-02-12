# Segmental Duplication Assembler (SDA)


# Download: #
```
git clone --recurse-submodules git://github.com/mvollger/SDA.git
```

# Install: #
The requirements for SDA are taken care of by two custom conda environments (sda-python-2 and sda-python-3). In order to run SDA, you must already have anaconda 3 installed on your system and you must be able to create conda environments. 

Once that is done, create `env_conda.cfg` so that it adds conda to your path. The `CONDA_PATH` variable must be set so that it points at you conda installation. 
In addition to your conda path we recommend adding gcc to your path and making your environment sparse.
We suggest using `gcc=6.4.0` because that is what we have used when installing SDA. 
Here is an example of what your `env_conda.cfg` might look like.
```
#!/bin/bash

# commands to clean your env...

# commands for adding conda to your path...  
export CONDA_PATH=/your/local/anaconda/install # e.g. /net/eichler/vol2/home/mvollger/anaconda3
export PATH=$CONDA_PATH/bin:$PATH

# commands to load gcc as an example I have provided the command I use to load gcc 
module load gcc/6.4.0

```

Once `env_conda.cfg` has been created and the anaconda environment is in your path (e.g. `source env_conda.cfg`), the `Makefile` can be run with:
```
make
```
Several people have run into an error with readToSNVList if it is complied with gcc 8.x so please complie with gcc 6. 




# Run: #

## Required input files: ## 
```
ref.fasta # collapsed representation of a segmental duplication
reads.orig.bam # reads aligning to the collapsed duplication. Note: if you want quiver/arrow error correction to happen these reads must be aligned using PacBio’s version of Blasr to preserve quality values. 
sda.config.json
```
The file `sda.config.json` must have four values set:
`MINCOV` should be set to a depth value above the average depth of sequencing error, `MAXCOV` should be set to a value around the average read depth, `MINTOTAL` should be set to the minimum expected coverage for a collapsed duplication, and `project` should just be a string identifier for your project. 
Below is an example of what the file might look like:
```
{
	"MINCOV" : 27, # Minimum depth over a PSV for it to be considered
	"MAXCOV" : 54, # Maximum depth over a PSV for it to be considered
	"MINTOTAL" : 83, # Minimum total depth at a PSV position for it to be considered. 
	"project": "CHM1_V4" # a name for your project (no spaces).
}
```
If you want to run with oxford nanopore technolgies (ONT) data you should add the following to `sda.config.json` above the `MINCOV` line: 
```
"ont" : "True",
```


### Optional input: ### 
If you are running SDA on a human genome you can provide additional files that will be used to compare your SDA results against. 
```
ref.fasta.bed  # bed file containing the locations of the duplications in GRCh38.
duplications.fasta # fasta file with the duplicated sequence from GRCh38. The header of each paralog must be in the format  “>{chromosome}:{start position}-{end position}”
```
If you include these files many additional outfiles are created. For example, `canu.summary.txt` which contains a summary of where the SDA contigs aligned to the reference and how well the match the reference. You can run the test case with and without `duplications.fasta` to see the difference. 


## Running: ##

SDA is run by a snakemake script. While SDA should run correctly as long as the required files are properly in place, you may find some knowledge of snakemake useful in running SDA.  

If the github repo is added to your path just type `SDA`; if not, type `/path/to/git/repo/SDA`.


## Output: ##
There are three main output files:
```
canu.assemblies.fasta
miniasm.assemblies.fasta
wtdbg.assemblies.fasta
```
These are fasta files containing the sequences for each paralog as determined by SDA. The prefixes say which assembler was used to generate the results - in my experience canu has worked best.

For a visualization of the correlation clustering results see `CC/mi.cuts.gml.pdf`

## Test cases: ##
There is a test case created by the make file in `TestCases/SDAtest`. Please test SDA on this before testing your own data (Note quiver will not run on this test case). 
A second test case in `TestCases/SDAtest2` will run quiver/arrow without error. 



### Common Java Error: ###
Some users will find that while creating the file `CC/mi.cuts.gml.pdf` they get the following error:
```
undefined symbol: FT_Done_MM_Var 
```
If you get this error we recommend re-installing `openjdk` as this has resolved the issue for others. 
```
source env_python3.cfg 
conda uninstall openjdk
conda install openjdk
conda install -c bioconda canu=1.8
```




# Identifying Collapsed Duplications #

I have written a snakemake (`ProcessCollapsedAssembly.py`) for identifying collapsed duplications within a de novo assembly. 
It requires the user to have a working install of `RepeatMasker` and setup a config script called `env_RM.cfg`, but otherwise the dependencies are taken care of by `sda-python-3`. An exmaple of my config is shown below:
```
module load perl/5.14.2
module load RepeatMasker/3.3.0
```
This process can be started by executing using this script `ProcessCollapsedAssembly.snake.sh` once the required input is in place. 
If you do not have a sungrid engine you will have to modify `ProcessCollapsedAssembly.snake.sh` to work on your cluster. 


## Required input: ##
This process only has one required input and that is a config file called `denovo.setup.config.json` which must be placed in a directory called config (`config/denovo.setup.config.json`).  An example of the `denovo.setup.config.json` file is shown below:
```
{	
	"asm" : "NA19240.fasta", # The denovo assembly to examine.
	"reads" : "reads.fofn", # A file of file names (FOFN) containing all the reads used in denovo assembly.
	"project" : "NA19240", # A project identifier, can be anything (no spaces). 
    "read_files_per_job" : 2  # the number of read files to submit to blasr at once. I recommend less than 15GB of data per job.
}
```
If you are running the pipeline on non-human data/assembly you should specify the `"species" : "your_species"`. `"your_species"` will be passed to RepeatMasker, thus it must be a vialid database for RepeatMasker. 

Once again `"ont" : "True"` can be added for use with ONT data. 
Additionally `"pbmm2" : "Ture"` can be added to use pbmm2 instead of blasr for alignments. 

There is a test case for this snakemake but it is rather large, 20Gb of data. To download it execute the following:
```
make TestCases/GenomeTest/ref.fasta
```

### Please Cite ###
Vollger, M. R. et al. Long-read sequence and assembly of segmental duplications. Nature Methods (2018). doi:10.1038/s41592-018-0236-3


