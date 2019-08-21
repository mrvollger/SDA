#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $DIR/../env_sda.sh
# rerunning dir in case ../env_sda.sh reset it 
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
envname=$DIR/sda-python-2
source activate $envname > /dev/null
export PATH=$envname/bin:$PATH

