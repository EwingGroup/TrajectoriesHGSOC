#!/bin/bash

# To run this script, do 
# qsub convertVCFtoBEDPE.sh  <VCF> 
#
# IDS is the bedpe to run clusterSV on
#
#$ -N convertVcfToBedpe
#$ -R y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=8G
#$ -pe sharedmem 4
#$ -l h_rt=300:00:00

unset MODULEPATH
. /etc/profile.d/modules.sh

module load igmm/apps/R/4.1.0

VCF=$1

SHORT=`basename $VCF .vcf`
echo $SHORT 
 
Rscript scripts/convertVCFtoBEDPE.R ${SHORT}.vcf

sed -i 's/chr//g' ${SHORT}.bedpe

