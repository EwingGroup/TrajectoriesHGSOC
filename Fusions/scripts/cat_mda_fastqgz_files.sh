#!/bin/bash 

# to run: qsub -t 1-n cat_mda_fastqgz_files.sh 

#$ -S /bin/bash
#$ -cwd
#$ -N catfastq 
#$ -tc 10
#$ -l h_vmem=1G
#$ -l h_rt=4:00:00
#$ -j y

PDIR=/exports/igmm/eddie/HGS-OvarianCancerA-SGP-WGS/fusions
SAMPLELIST=$PDIR/MDA_ids.list
sampleID=`sed -n "${SGE_TASK_ID}p" $SAMPLELIST`

DATDIR=/exports/igmm/eddie/HGS-OvarianCancerA-SGP-WGS/EGA_MDA/rnaseq
OUTDIR=$PDIR/MDA_fastq_merged
files1=`ls $DATDIR/"$sampleID"-*_R1.fastq.gz`
echo "files1 is $files1"
cat $files1 >  $OUTDIR/"$sampleID"_R1.fastq.gz
files2=`ls $DATDIR/"$sampleID"-*_R2.fastq.gz`
echo "files2 is $files2"
cat $files2 >  $OUTDIR/"$sampleID"_R2.fastq.gz
