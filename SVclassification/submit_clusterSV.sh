#!/bin/bash

# To run this script, do 
# qsub submit_clusterSV.sh  <IDS> 
#
# IDS is the bedpe to run clusterSV on
#
#$ -N clusterSV
#$ -R y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=16G
#$ -l h_rt=24:00:00

unset MODULEPATH
. /etc/profile.d/modules.sh

module load igmm/apps/R/4.1.0
IDS=$1

GZVCF=`head -n $SGE_TASK_ID $IDS | tail -n 1 | awk '{ print $1 }'`
SHORT=`basename $GZVCF .vcf`
echo $SHORT 
 

Rscript ../scripts/convertVCFtoBEDPE.R ${SHORT}.vcf

sed -i 's/chr//g' ${SHORT}.bedpe

Rscript ../scripts/clustering_index.R ${SHORT}.bedpe 16

