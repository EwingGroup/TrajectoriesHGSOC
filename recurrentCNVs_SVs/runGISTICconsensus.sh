#!/bin/bash

# To run this script, do 
# qsub -t 1-n submit_example_batcharray_inR.sh <IDS>
#
# IDS is the list of ids to run the R script on
#
#$ -N gistic_consensus
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=8G
#$ -l h_rt=24:00:00

unset MODULEPATH
. /etc/profile.d/modules.sh

module load igmm/apps/gistic/2.0.23

./gistic2 \
-b HGSOC_consensus \
-refgene /exports/igmm/eddie/bioinfsvice/ggrimes/gistic/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat \
-seg consensus.seg \
-ta 0.1 \
-td 0.1 \
-qvt 0.05 \
-conf 0.95 \
-fname hgsoc \
-broad 1 \
-maxseg 4000 
