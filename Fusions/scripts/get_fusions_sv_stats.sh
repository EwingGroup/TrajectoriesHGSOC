#!/bin/bash 

# run this via: qsub -v COHORT=TCGA_results -v IDLIST=TCGA_ids.list -v FUSIONLIST=fusions.list -v OUTPUT=TCGA_fusions_filt.txt get_fusions_sv_stats.sh 
# can also run on command line: 
# get_fusions_sv_stats.sh COHORT IDLIST FUSIONLIST OUTPUT

#$ -S /bin/bash
#$ -cwd
#$ -N filtfusions
#$ -l h_vmem=1G
#$ -l h_rt=01:00:00
#$ -j y

if [ -z $COHORT ]; then  
	COHORT=$1 # TCGA_results, etc
	IDLIST=$2
	FUSIONLIST=$3 # This is a list of FusionNames that we are interested in (ie SAA4--SAA1) 
	OUTPUT=$4 # The output is a table with a row for each FusionName and columns that have the number of patients in the cohort that have the fusion detected at various filtering stages 
fi 
echo COHORT is $COHORT, IDLIST is $IDLIST, FUSIONLIST is $FUSIONLIST, and OUTPUT is $OUTPUT
myfilters=(ffpm mantasv mantap lumpysv lumpyp)

echo "FusionName filtered.FFPM MantaSV MantaSVpass LumpySV LumpySVpass " | tr ' ' '\t' > $OUTPUT
for fusion in `cat $FUSIONLIST`; do
	outline=$fusion
	declare -A filt_array
	for myfilt in "${myfilters[@]}"; do 
		filt_array[$myfilt]=0
	done 
	for id in `cat $IDLIST`; do
		myfile=`ls $COHORT/$id/star-fusion.fusion_candidates.preliminary.filtered.FFPM.SVcalls`
		x1=`grep $fusion $myfile | wc -l`
		if [ $x1 -gt 0 ]; then
			(( filt_array[ffpm]++ ))
		
			x=`grep $fusion $myfile | cut -f6-9 | grep -v "0:0::.	0:0::.	0:0::.	0:0::." | wc -l`
			if [ $x -gt 0 ]; then 
				(( filt_array[mantasv]++ ))
			fi 
			x=`grep $fusion $myfile | cut -f6-9 | grep PASS | wc -l`
			if [ $x -gt 0 ]; then 
				(( filt_array[mantap]++ ))
			fi 
			x=`grep $fusion $myfile | cut -f10-13 | grep -v "0:0::.	0:0::.	0:0::.	0:0::." | wc -l`
			if [ $x -gt 0 ]; then 
				(( filt_array[lumpysv]++ ))
			fi 
			x=`grep $fusion $myfile | cut -f10-13 | grep PASS | wc -l`
			if [ $x -gt 0 ]; then 
				(( filt_array[lumpyp]++ ))
			fi
		fi
	done
	for myfilt in "${myfilters[@]}"; do 
		outline=$outline" "${filt_array[$myfilt]}
	done
	echo $outline | tr ' ' '\t' >> $OUTPUT
done



