#!/bin/bash 

#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=10G 
#$ -j y

. /etc/profile.d/modules.sh
MODULEPATH=$MODULEPATH:/exports/igmm/software/etc/el7/modules 

module load anaconda/5.3.1
module load igmm/apps/nextflow/20.12.0-edge
module load igmm/apps/hmmer/3.1b2

nextflow run nf-core/rnafusion/build-ctat.nf \
	--outdir references \
	--fasta /exports/igmm/eddie/NextGenResources/nextflow/GRCh38_files/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	--gtf /exports/igmm/eddie/NextGenResources/nextflow/GRCh38_files/Homo_sapiens.GRCh38.103.gtf \
	-profile conda \
	-c conf/eddie2.config
