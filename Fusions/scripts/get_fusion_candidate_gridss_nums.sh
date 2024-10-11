#!/bin/bash

COHORT=$1 #TCGA_results, etc
IDLIST=$2
OUTPUT=$3
echo "sampleID	filtered.FFPM	gridssSV	gridSVPass" > $OUTPUT

for id in `cat $IDLIST`; do
	outline=$id
	myfile=`ls $COHORT/$id/star-fusion.fusion_candidates.preliminary.filtered.FFPM.gridssSVcalls` 
	x=`sed 1d $myfile | wc -l`
	grid=`sed 1d $myfile | cut -f6-9 | grep -v "0:0::.	0:0::.	0:0::.	0:0::." | wc -l`
	gridpass=`sed 1d $myfile | cut -f6-9 | grep PASS | wc -l`
	outline="$outline $x $grid $gridpass"
	echo $outline >> $OUTPUT
done
cat $OUTPUT | tr ' ' '\t' > tmp
mv tmp $OUTPUT

