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
IDLIST=$PROJDIR/fusions/MDA_wgs_ids.list
FUSION_DIR=$PROJDIR/fusions/starfusion_results/MDA_results
SV_DIR=$PROJDIR/variants/bcbio

scratch=~/scratch/mergesv/$JOB_ID.$SGE_TASK_ID
mkdir -p $scratch
sampleID=`sed -n "${SGE_TASK_ID}p" $IDLIST | cut -f1`

fusionbn=star-fusion.fusion_candidates.preliminary.filtered.FFPM
fusionfile=$FUSION_DIR/$sampleID/star-fusion.preliminary/$fusionbn
OUTDIR=$FUSION_DIR/$sampleID

varid=`sed -n "${SGE_TASK_ID}p" $IDLIST | cut -f2`
for svcaller in manta; do
	SVfile=$SV_DIR/$varid/*/$varid"-$svcaller.vcf.gz"
	python $PROJDIR/fusions/scripts/search_for_sv_breakpoint_in_fusions.py $fusionfile $SVfile $scratch/$fusionbn.$svcaller
done
echo "FusionName JunctionReadCount SpanningFragCount LeftBreakpoint RightBreakpoint MantaSV_Left MantaSV_Right" | tr ' ' '\t' > $OUTDIR/$fusionbn.SVcalls
cat $scratch/$fusionbn.manta >> $OUTDIR/$fusionbn.SVcalls

