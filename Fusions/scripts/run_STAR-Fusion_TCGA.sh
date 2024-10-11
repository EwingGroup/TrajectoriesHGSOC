#!/bin/bash 

# run this via qsub -t 1-n run_STAR-Fusion_TCGA.sh 

#$ -S /bin/bash
#$ -cwd
#$ -N starfusion
#$ -tc 10
#$ -l h_vmem=6G
#$ -l h_rt=24:00:00
#$ -pe sharedmem 8
#$ -j y

. /etc/profile.d/modules.sh 
MODULEPATH=$MODULEPATH:/exports/igmm/software/etc/el7/modules 
module load igmm/apps/perl/5.24.0
module load igmm/apps/STARFusion/1.10.0
module load igmm/apps/samtools/1.6
module load igmm/apps/STAR/2.7.8a 
## Add STAR version 2.7.8a to your path 
# PATH=$PATH:/gpfs/igmmfs01/eddie/NextGenResources/software/STAR_2.7.8a/STAR-2.7.8a/bin/Linux_x86_64

CTAT_GENOME_LIB=/exports/igmm/eddie/HGS-OvarianCancerA-SGP-WGS/fusions/STAR-Fusion/GRCh38_gencode_v37_CTAT_lib_Mar012021.source/ctat_genome_lib_build_dir
# CTAT_GENOME_LIB=/exports/igmm/software/pkg/el7/apps/STARFusion/1.10.0/ctat-genome-lib-builder/GRCh38_gencode_v37_CTAT_lib_Mar012021.source/ctat_genome_lib_build_dir
# starbin=/exports/igmm/eddie/HGS-OvarianCancerA-SGP-WGS/fusions/STAR-Fusion/STAR-Fusion-v1.10.0
datadir=/exports/igmm/eddie/HGS-OvarianCancerA-SGP-WGS/TCGA_US_OV/rnaseq/reads
SAMPLELIST=/exports/igmm/eddie/HGS-OvarianCancerA-SGP-WGS/fusions/TCGA_ids.list

sampleID=`sed -n "${SGE_TASK_ID}p" $SAMPLELIST`

echo "sampleID is $sampleID"

STAR-Fusion --genome_lib_dir $CTAT_GENOME_LIB \
	--left_fq $datadir/$sampleID"_R1.fastq.gz" \
	--right_fq $datadir/$sampleID"_R2.fastq.gz" \
	--output_dir TCGA_results/$sampleID

