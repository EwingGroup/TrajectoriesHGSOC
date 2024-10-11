#!/bin/bash 

# to run: qsub -t 1-n cat_AOCS_fastqgz_files.sh 

#$ -S /bin/bash
#$ -cwd
#$ -N catfastq 
#$ -tc 10
#$ -l h_vmem=2G
#$ -l h_rt=6:00:00
#$ -j y

PDIR=/exports/igmm/eddie/HGS-OvarianCancerA-SGP-WGS/fusions
SAMPLELIST=$PDIR/AOCS_ids_tomerge.list
sampleID=`sed -n "${SGE_TASK_ID}p" $SAMPLELIST | awk '{print $1}'`

DATDIR=/exports/igmm/eddie/HGS-OvarianCancerA-SGP-WGS/rnaseq/raw_reads/AOCS
OUTDIR=$PDIR/AOCS_fastq_merged
files1=`ls $DATDIR/"$sampleID"_*_R1.fastq.gz`
echo "files1 is $files1"
echo "cat $files1 >  $OUTDIR/"$sampleID"_R1.fastq.gz"
cat $files1 >  $OUTDIR/"$sampleID"_R1.fastq.gz
nlines=`gzip -dc $OUTDIR/"$sampleID"_R1.fastq.gz | wc -l`
echo "number of lines in R1 is $nlines" 
files2=`ls $DATDIR/"$sampleID"_*_R2.fastq.gz`
echo "files2 is $files2"
cat $files2 >  $OUTDIR/"$sampleID"_R2.fastq.gz
nlines=`gzip -dc $OUTDIR/"$sampleID"_R2.fastq.gz | wc -l`
echo "number of lines in R2 is $nlines" 

