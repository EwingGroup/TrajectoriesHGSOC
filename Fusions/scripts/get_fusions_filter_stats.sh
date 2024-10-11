#!/bin/bash 
# run this via: qsub -v COHORT=TCGA_results -v IDLIST=TCGA_ids.list -v FUSIONLIST=fusions.list -v OUTPUT=TCGA_fusions_filtcnt.txt get_fusions_filter_stats.sh 
# can also run on command line: 
# get_fusions_filter_stats.sh COHORT IDLIST FUSIONLIST OUTPUT

#$ -S /bin/bash
#$ -cwd
#$ -N filtfusions
#$ -l h_vmem=1G
#$ -l h_rt=40:00:00
#$ -j y

if [ -z $COHORT ]; then 
	echo "what?!"
	COHORT=$1 # TCGA_results, etc
	IDLIST=$2
	FUSIONLIST=$3 # This is a list of FusionNames that we are interested in (ie SAA4--SAA1) 
	OUTPUT=$4 # The output is a table with a row for each FusionName and columns that have the number of patients in the cohort that have the fusion detected at various filtering stages 
fi 

echo COHORT is $COHORT, IDLIST is $IDLIST, FUSIONLIST is $FUSIONLIST, and OUTPUT is $OUTPUT
mysuffixes=(preliminary filtered.FFPM   annot_filter.pass   RTartifact.pass minFFPM.0.1.pass)
echo "FusionName preliminary filtered.FFPM annot_filter.pass RTartifact.pass minFFPM.0.1.pass" | tr ' ' '\t' > $OUTPUT
for fusion in `cat $FUSIONLIST`; do
	outline=$fusion 
	declare -A filt_array
	for mysuf in "${mysuffixes[@]}"; do
		filt_array[$mysuf]=0
		for id in `cat $IDLIST`; do 
			myfile=`ls $COHORT/$id/star-fusion.preliminary/star-fusion.fusion_candidates.*$mysuf`
			x=`grep $fusion $myfile | wc -l`
			if [ $x -gt 0 ]; then 
				(( filt_array[$mysuf]++ ))
			fi
		done
		outline=$outline" "${filt_array[$mysuf]}
	done
	echo $outline | tr ' ' '\t' >> $OUTPUT
done



