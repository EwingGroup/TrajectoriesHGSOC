for i in *consensus.vcf
do

sample=`grep '#CHROM' $i | cut -f 11`
grep -v '#' $i > ${sample}.vars_orig
sed "s/$/\t${sample}/" ${sample}.vars_orig > new.${sample}.vars

done

Rscript scripts/convertVCFtoBEDPE.R ${SHORT}.vcf

sed -i 's/chr//g' ${SHORT}.bedpe