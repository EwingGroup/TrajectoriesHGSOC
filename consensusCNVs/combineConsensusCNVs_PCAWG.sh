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
PCAWG=`head -n $SGE_TASK_ID $PARAMS | tail -n 1 | awk '{ print $2 }'`
OUTDIR=$2
OVERLAP=$3

CONSENSUS="/exports/igmm/eddie/HGS-OvarianCancerA-SGP-WGS/SVhotspots/consensus_cnvs/${SAMPLE}-cnvkit_CLImaT_${OVERLAP}_all_lengths-uniq.bed"

#Lose PCAWG header
grep -v 'chromosome' $PCAWG | cut -f 1,2,3,4 > $OUTDIR/${SAMPLE}-pcawg.bed 

bedtools intersect -a ${CONSENSUS} -b ${OUTDIR}/${SAMPLE}-pcawg.bed -wo -f $OVERLAP -r > $OUTDIR/${SAMPLE}-consensus_pcawg_tmp.bed

Rscript --no-save scripts/combineCNVkit_CLImAT.R ${OUTDIR}/${SAMPLE}-consensus_pcawg_tmp.bed ${OUTDIR}/${SAMPLE}-consensus_pcawg_${OVERLAP}_all_lengths-uniq.bed



rm ${OUTDIR}/${SAMPLE}-consensus_pcawg_tmp.bed ${OUTDIR}/${SAMPLE}-pcawg.bed




