#!/bin/bash 

# to run: qsub -t 1-n cat_fastqgz_files.sh 

#$ -S /bin/bash
#$ -cwd
#$ -N catfastq 
#$ -tc 10
#$ -l h_vmem=1G
#$ -l h_rt=4:00:00
#$ -j y

PDIR=/exports/igmm/eddie/HGS-OvarianCancerA-SGP-WGS/rnaseq/raw_reads/SHGSOC
SAMPLELIST=$PDIR/renamed_merged/usable_ids.txt
sampleID=`sed -n "${SGE_TASK_ID}p" $SAMPLELIST`

DATDIR=$PDIR/renamed
OUTDIR=$PDIR/renamed_merged
files1=`ls $DATDIR/"$sampleID"_*R1_001.fastq.gz`
echo "files1 is $files1"
cat $files1 >  $OUTDIR/"$sampleID"_R1.fastq.gz
files2=`ls $DATDIR/"$sampleID"_*R2_001.fastq.gz`
cat $files2 >  $OUTDIR/"$sampleID"_R2.fastq.gz
