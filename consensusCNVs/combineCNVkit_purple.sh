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
PURPLE=`head -n $SGE_TASK_ID $PARAMS | tail -n 1 | awk '{ print $4 }'`
OUTDIR=$2
OVERLAP=$3


tail -n +14 $CLIMAT | cut -f 1,2,3,4,5,7 | awk '$4 !=2' | awk '$6>50'> $OUTDIR/${SAMPLE}-CLImAT.bed  

#Intersect with CNVkit
grep -v 'chrUn' $CNVKIT |grep -v 'chrY'|grep -v 'random'|grep -v 'M'|sed 's/chr//g'|sed 's/X/23/g'|awk 'BEGIN { OFS = "\t" } { $7 = $3 - $2 } 1' > $OUTDIR/${SAMPLE}-cnvkit-nopseudo_all_lengths.bed

#Format Purple
grep -v 'chromosome' $PURPLE |sed 's/chr//g'  | awk '$10!="STRUCTURAL_VARIANT"'|sed 's/X/23/g' >  $OUTDIR/${SAMPLE}-purple.bed



bedtools intersect -a ${OUTDIR}/${SAMPLE}-cnvkit-nopseudo_all_lengths.bed -b ${OUTDIR}/${SAMPLE}-purple.bed -wo -f $OVERLAP  > $OUTDIR/${SAMPLE}-cnvkit-purple_tmp.bed

Rscript --no-save scripts/combineCNVkit_CLImAT.R ${OUTDIR}/${SAMPLE}-cnvkit-purple_tmp.bed ${OUTDIR}/${SAMPLE}-cnvkit-purple_${OVERLAP}_all_lengths.bed


rm ${OUTDIR}/${SAMPLE}-cnvkit-purple_tmp.bed  




