#!/bin/bash

# To run this script, do 
# qsub -t 1-n submit_example_batcharray_inR.sh <IDS>
#
# IDS is the list of ids to run the R script on
#
#$ -N gistic
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -pe sharedmem 16
#$ -l h_vmem=4G
#$ -l h_rt=24:00:00

unset MODULEPATH
. /etc/profile.d/modules.sh

module load igmm/apps/gistic/2.0.23

./gistic2 \
-b HGSOC_purple_x \
-refgene /exports/igmm/eddie/bioinfsvice/ggrimes/gistic/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat \
-seg all.seg.purple \
-ta 0.1 \
-td 0.1 \
-qvt 0.05 \
-conf 0.95 \
-maxseg 4000 \
-fname hgsoc_purple_x \
-broad 1 \
-savegene 1 \
-rx 0
