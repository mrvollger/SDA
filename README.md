# Segmental Duplication Assembler (SDA)


# Download: #
```
git clone --recurse-submodules git://github.com/mvollger/SDA.git
```

# Install: #
The requirements for SDA are taken care of by two custom conda environments (sda-python-2 and sda-python-3). In order to run SDA, you must already have anaconda 3 installed on your system and you must be able to create conda environments. 

Once that is done, modify `env_conda.cfg` so that it adds conda to your path. The `CONDA_PATH` variable must be set so that it points at you conda installation. Here is an example of what your `env_conda.cfg` might look like.
```
#!/bin/bash

export CONDA_PATH=/your/local/anaconda/install # e.g. /net/eichler/vol2/home/mvollger/anaconda3
export PATH=$CONDA_PATH/bin:$PATH

```

Once `env_conda.cfg` has been updated and the anaconda environment is in your path (e.g. `source env_conda.cfg`), the `Makefile` can be run with:
```
make
```
We suggest using `gcc=7.x` because that is what we have used when installing SDA. 




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

## Test case: ##
There is a test case created by the make file in `TestCases/SDAtest`. Please test SDA on this before testing your own data (Note quiver will not run on this test case). 



### Common Java Error: ###
Some users will find that while creating the file `CC/mi.cuts.gml.pdf` they get the following error:
```
undefined symbol: FT_Done_MM_Var 
```
If you get this error we recommend re-installing the `java-jdk` as this has resolved the issue for others. 
```
source activate sda-python-3 
conda uninstall java-jdk
conda install java-jdk
```




# Identifying Collapsed Duplications #

I have written a snakemake (`ProcessCollapsedAssembly.py`) for identifying collapsed duplications within a de novo assembly. It requires the user to have a working install of `RepeatMasker`, but otherwise the dependencies are taken care of by `sda-python-3`. 

This process can be started by executing using this script `ProcessCollapsedAssembly.snake.sh` once the required input is in place. 

## Required input: ##
This process only has one required input and that is a config file called `denovo.setup.config.json` which must be placed in a directory called config (`config/denovo.setup.config.json`).  An example of the `denovo.setup.config.json` file is shown below:
```
{	
	"asm" : "NA19240.fasta", # The denovo assembly to examine.
	"reads" : "reads.fofn", # A file of file names (FOFN) containing all the reads used in denovo assembly.
	"reference " :  "~mvollger/assemblies/hg38/ucsc.hg38.no_alts.fasta", # a path to a local download of UCSC’s hg38. 
        "genes" : "/net/eichler/vol2/home/mvollger/assemblies/hg38/hg38.gene.locations.bed",
	"project" : "NA19240", # A project identifier, can be anything (no spaces). 
    	"bax_per_job" : 10  # the number of read files to submit to blasr at once. I recommend less than 15. 
}
```
Once again `"ont" : "True"` can be added for use with ONT data. 


