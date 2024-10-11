#!/bin/bash

# To run this script, do 
# qsub -t 1-324 getConsensusCalls.sh  getConsensusCalls.py <params_consensus_calls> 
#
#
#$ -N consensusSV
#$ -R y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=8G
#$ -l h_rt=01:00:00

unset MODULEPATH
. /etc/profile.d/modules.sh

module load anaconda

IDS=$1

manta=`head -n $SGE_TASK_ID $IDS | tail -n 1 | awk '{ print $1 }'`
gridss=`head -n $SGE_TASK_ID $IDS | tail -n 1 | awk '{ print $2 }'`
sample=`head -n $SGE_TASK_ID $IDS | tail -n 1 | awk '{ print $3 }'`

source activate viola-sv

python ../scripts/getConsensusCalls.py $manta $gridss $sample

source deactivate
