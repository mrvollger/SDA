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
We suggest using `gcc=7.x` because that is what we have used we installing SDA. 




# Run: #

## Required input files: ## 
```
ref.fasta # collapsed representation of a segmental duplication
reads.orig.bam # reads aligning to the collapsed duplication. Note: if you want quiver/arrow error correction to happen these reads must be aligned using PacBio’s version of Blasr to preserve quality values. 
sda.config.json
```
The file `sda.config.json` must have three values set:
 `MINCOV` should be set to a depth value above the average depth of sequencing error, `MAXCOV` should be set to a value around the average read depth, and `MINTOTAL` should be set to the minimum expected coverage for a collapsed duplication. 
Below is an example of what the file might look like:
```
{
	"MINCOV" : 27, # Minimum depth over a PSV for it to be considered
	"MAXCOV" : 54, # Maximum depth over a PSV for it to be considered
	"MINTOTAL" : 83, # Minimum total depth at a PSV position for it to be considered. 
}
```



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



