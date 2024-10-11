#!/bin/bash 

#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=10G 
#$ -j y

. /etc/profile.d/modules.sh
MODULEPATH=$MODULEPATH:/exports/igmm/software/etc/el7/modules 

module load anaconda/5.3.1
module load igmm/apps/nextflow/20.04.1

nextflow run my-build-ctat.nf \
	--outdir references \
	--fasta /exports/igmm/eddie/NextGenResources/nextflow/GRCh38_files/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	--gtf /exports/igmm/eddie/NextGenResources/nextflow/GRCh38_files/Homo_sapiens.GRCh38.103.gtf \
	-profile conda \
	-c conf/eddie2.config
