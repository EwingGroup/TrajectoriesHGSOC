#!/bin/bash 

# Run this via: qsub -t1-n run_merge_sv_with_fusions.sh 

#$ -S /bin/bash
#$ -cwd
#$ -N mergewSV
#$ -tc 10
#$ -l h_vmem=1G
#$ -l h_rt=00:30:00
#$ -j y

. /etc/profile.d/modules.sh 
MODULEPATH=$MODULEPATH:/exports/igmm/software/etc/el7/modules
module load roslin/python/3.6.8  # This version of python has pyvcf installed

PROJDIR=/exports/igmm/eddie/HGS-OvarianCancerA-SGP-WGS
IDLIST=$PROJDIR/fusions/SHGSOC_usable_ids.list
FUSION_DIR=$PROJDIR/fusions/SHGSOC_results
SV_DIR=$PROJDIR/variants/bcbio

scratch=~/scratch/mergesv/$JOB_ID.$SGE_TASK_ID
mkdir -p $scratch
sampleID=`sed -n "${SGE_TASK_ID}p" $IDLIST`

fusionbn=star-fusion.fusion_candidates.preliminary.filtered.FFPM
fusionfile=$FUSION_DIR/$sampleID/star-fusion.preliminary/$fusionbn
OUTDIR=$FUSION_DIR/$sampleID

if [ `echo $sampleID | grep R2` ]; then 
	varid=`echo $sampleID | sed 's/R2//'`
else
	varid=`echo $sampleID | sed 's/R//'`
fi 

for svcaller in manta lumpy; do
	SVfile=$SV_DIR/$varid/$varid"T"/$varid-$svcaller.vcf.gz
	python $PROJDIR/fusions/scripts/search_for_sv_breakpoint_in_fusions.py $fusionfile $SVfile $scratch/$fusionbn.$svcaller
	SVfile=$SV_DIR/$varid/$varid"N"/$varid"N"-$svcaller.vcf.gz
	python $PROJDIR/fusions/scripts/search_for_sv_breakpoint_in_fusions.py $fusionfile $SVfile $scratch/$fusionbn.$svcaller"-N"
done
echo "FusionName JunctionReadCount SpanningFragCount LeftBreakpoint RightBreakpoint MantaSV_Left MantaSV_Right MantaNSV_Left MantaNSV_Right LumpySV_Left LumpySV_Right LumpyNSV_Left LumpyNSV_Right" | tr ' ' '\t' > $OUTDIR/$fusionbn.SVcalls

paste $scratch/$fusionbn.manta <(cut -f6- $scratch/$fusionbn.manta"-N") <(cut -f6- $scratch/$fusionbn.lumpy) <(cut -f6- $scratch/$fusionbn.lumpy"-N") >> $OUTDIR/$fusionbn.SVcalls 

