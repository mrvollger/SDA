#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# loads conda enviormaent and the CONDA_PATH
source $DIR/env_conda.sh

envname=$DIR/envs/sda-python-2
source activate $envname > /dev/null
export PATH=$envname/bin:$PATH

