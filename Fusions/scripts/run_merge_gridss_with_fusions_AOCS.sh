#!/bin/bash 

# Run this via: qsub -t1-n run_merge_gridss_with_fusions.sh 

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
module load igmm/libs/htslib/1.13 # for tabix

PROJDIR=/exports/igmm/eddie/HGS-OvarianCancerA-SGP-WGS
SDIR=$PROJDIR/fusions/scripts
IDLIST=$PROJDIR/fusions/AOCS_usable_ids.list
FUSION_DIR=$PROJDIR/fusions/AOCS_results
SV_DIR=$PROJDIR/variants/gridss-purple-linx/variants

scratch=~/scratch/mergesv/$JOB_ID.$SGE_TASK_ID
mkdir -p $scratch
sampleID=`sed -n "${SGE_TASK_ID}p" $IDLIST`

fusionbn=star-fusion.fusion_candidates.preliminary.filtered.FFPM
fusionfile=$FUSION_DIR/$sampleID/star-fusion.preliminary/$fusionbn
OUTDIR=$FUSION_DIR/$sampleID

varid=$sampleID"_PrimaryTumour" 
SVfile=$SV_DIR/$varid/gridss/$varid.gridss.vcf.gz
tabix -f $SVfile
python $SDIR/search_for_sv_breakpoint_in_fusions.py $fusionfile $SVfile $scratch/$fusionbn.gridss
SVfile=$SV_DIR/$varid/gridss/$varid.gridss.full.somatic.vcf.gz
gzip -dc $SVfile | grep -v ^"##ALT" > $scratch/somatic.vcf
bgzip $scratch/somatic.vcf
tabix -f $scratch/somatic.vcf.gz
python $SDIR/search_for_sv_breakpoint_in_fusions.py $fusionfile $scratch/somatic.vcf.gz $scratch/$fusionbn.gridSomatic

echo "FusionName JunctionReadCount SpanningFragCount LeftBreakpoint RightBreakpoint Gridss_Left Gridss_Right Gridsom_Left Gridsom_Right" | tr ' ' '\t' > $OUTDIR/$fusionbn.gridssSVcalls

paste $scratch/$fusionbn.gridss <(cut -f6- $scratch/$fusionbn.gridSomatic) >> $OUTDIR/$fusionbn.gridssSVcalls 

