#!/usr/env/bin bash
module purge
. /etc/profile.d/modules.sh
module load modules modules-init modules-gs/prod modules-eichler
module load anaconda/201710
# properties = {properties}
{exec_job}
