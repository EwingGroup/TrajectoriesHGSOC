#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l h_rt=04:00:00
#$ -l h_vmem=4G
#$ -N combineCNVkit_CLImAT

. /etc/profile.d/modules.sh
MODULEPATH=$MODULEPATH:/exports/igmm/software/etc/el7/modules

module load igmm/apps/R/3.6.0 
module load igmm/apps/BEDTools/2.27.1

PARAMS=$1
SAMPLE=`head -n $SGE_TASK_ID $PARAMS | tail -n 1 | awk '{ print $1 }'`
CLIMAT=`head -n $SGE_TASK_ID $PARAMS | tail -n 1 | awk '{ print $2 }'`
CNVKIT=`head -n $SGE_TASK_ID $PARAMS | tail -n 1 | awk '{ print $3 }'`
PCAWG=`head -n $SGE_TASK_ID $PARAMS | tail -n 1 | awk '{ print $5 }'`
OUTDIR=$2
OVERLAP=$3


tail -n +14 $CLIMAT | cut -f 1,2,3,4,5,7 | awk '$4 !=2' | awk '$6>50'> $OUTDIR/${SAMPLE}-CLImAT.bed  

#Lose PCAWG header
grep -v 'chromosome' $PCAWG | cut -f 1,2,3,4 > $OUTDIR/${SAMPLE}-pcawg.bed 

bedtools intersect -a ${OUTDIR}/${SAMPLE}-CLImAT.bed -b ${OUTDIR}/${SAMPLE}-pcawg.bed -wo -f $OVERLAP  > $OUTDIR/${SAMPLE}-climat-pcawg.bed

Rscript --no-save scripts/combineClimat_PCAWG.R ${OUTDIR}/${SAMPLE}-climat-pcawg.bed ${OUTDIR}/${SAMPLE}-climat-pcawg_${OVERLAP}_all_lengths.bed


#rm ${OUTDIR}/${SAMPLE}-climat-pcawg.bed  




