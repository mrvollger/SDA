# Segmental Duplication Assembler (SDA)


# Download: #
```
git clone --recurse-submodules -j8 git://github.com/mvollger/SDA.git
```

# Install: #
The requirements for SDA are taken care of my two conda environments (abp-python-2 and abp-python-3). In order to get this running you must already have anaconda 3 installed somewhere on your system. 

Once that is done modify `env_conda.cfg` so that it adds conda to your path. The `CONDA_PATH` variable must be set so that it points at you conda installation. Here is an example of what your `env_conda.cfg` might look like.
```
#!/bin/bash
unset PYTHONPATH

# unsetting your current env can be important, this is how I unload mine, but it likely does not translate to your system
module purge
. /etc/profile.d/modules.sh
module load modules modules-init modules-gs/prod modules-eichler

# conda path should be updated by the end user
# MUST BE CHANGED BY USER!!!
export CONDA_PATH=/your/local/anaconda/install # e.g. /net/eichler/vol2/home/mvollger/anaconda3
export PATH=$CONDA_PATH/bin:$PATH
```

Once the anaconda environment is available the `Makefile` can be run with:
```
make
```

# Setup: #
SDA requires some files to be run.
```
ref.fasa # collapsed representation of a segmental duplication
reads.orig.bam # reads aligning to the collapsed duplication. Note if you want error correction by quiver to happen these reads must be aligned using PacBio’s version of Blasr to preserve quality values. 

```


# Run: #

SDA is run by a snakemake script. While SDA should run correctly as long as the setup files are properly in place you may find some knowledge of snakemake useful in running SDA.  

If the github repo is added to your path just type `SDA` if not type `/path/to/git/repo/SDA`.



