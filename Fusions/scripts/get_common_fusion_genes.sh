#!/bin/bash

COHORT=$1
IDLIST=$2
OUTPUT=$3
echo -e "SampleID\tFusionName\tLeftGene\tRightGene\tannots" > $OUTPUT
for id in `cat $IDLIST`; do 
	myfile=$COHORT/$id/star-fusion.fusion_predictions.abridged.tsv
	cut -f1,7,9,17 $myfile | grep -v LeftGene | sort -u | awk '{print id"\t"$1"\t"$2"\t"$3"\t"$4}' id=$id >> $OUTPUT
done

