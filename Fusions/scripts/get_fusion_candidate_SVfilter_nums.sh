#!/bin/bash

COHORT=$1 #TCGA_results, etc
IDLIST=$2
OUTPUT=$3
echo "sampleID	filtered.FFPM	MantaSV	MantaSVpass LumpySV LumpySVpass GridssSV	GridssSVpass" > $OUTPUT

for id in `cat $IDLIST`; do
	outline=$id
	myfile=$COHORT/$id/star-fusion.fusion_candidates.preliminary.filtered.FFPM.SVcalls 
	if [ -f "$myfile" ]; then 
		x=`sed 1d $myfile | wc -l`
		mnta=`sed 1d $myfile | cut -f6-9 | grep -v "0:0::.	0:0::.	0:0::.	0:0::." | wc -l`
		mntapass=`sed 1d $myfile | cut -f6-9 | grep PASS | wc -l`
		lpy=`sed 1d $myfile | cut -f10-13 | grep -v "0:0::.	0:0::.	0:0::.	0:0::." | wc -l`
		lpypass=`sed 1d $myfile | cut -f10-13 | grep PASS | wc -l`
	fi
	myfile=$COHORT/$id/star-fusion.fusion_candidates.preliminary.filtered.FFPM.gridssSVcalls
	if [ -f "$myfile" ]; then 
		grid=`sed 1d $myfile | cut -f6-9 | grep -v "0:0::.  0:0::.  0:0::.  0:0::." | wc -l`
    	gridpass=`sed 1d $myfile | cut -f6-9 | grep PASS | wc -l`
	else
		grid="NA"
		gridpass="NA"
	fi 
	outline="$outline $x $mnta $mntapass $lpy $lpypass $grid $gridpass"
	echo $outline >> $OUTPUT
done
cat $OUTPUT | tr ' ' '\t' > tmp
mv tmp $OUTPUT
