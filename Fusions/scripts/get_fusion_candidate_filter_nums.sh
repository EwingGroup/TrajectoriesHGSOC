#!/bin/bash

COHORT=$1 #TCGA_results, etc
IDLIST=$2
OUTPUT=$3
mysuffixes=(preliminary	filtered.FFPM	annot_filter.pass	RTartifact.pass	minFFPM.0.1.pass)
echo "sampleID	preliminary	filtered.FFPM	annot_filter.pass	RTartifact.pass	minFFPM.0.1.pass" > $OUTPUT

for id in `cat $IDLIST`; do
	outline=$id
	for mysuf in "${mysuffixes[@]}"; do 
		myfile=`ls $COHORT/$id/star-fusion.preliminary/star-fusion.fusion_candidates.*$mysuf`
		x=`wc -l < $myfile`
		outline=$outline" "$x
	done
	echo $outline >> $OUTPUT
done
cat $OUTPUT | tr ' ' '\t' > tmp
mv tmp $OUTPUT
