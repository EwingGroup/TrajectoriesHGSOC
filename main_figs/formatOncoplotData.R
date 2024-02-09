##################################
# Combine all data for oncoplots #
##################################

library(tidyverse)
library(data.table)


#############################
# Format SNV and INDEL data #
#############################
snvsParsed <- fread("data/hgsocSNVs.txt.gz")

snvs <- snvsParsed %>%
  
  # format sample names, format SIFT and PolyPhen predictions and annotate pathogenic variants
  mutate(V11 = gsub("-ensemble-annotated-oxog_filtered.vep107.vcf.gz", "", V11),
         V11 = gsub("_2$|-P$", "", V11),
         V5 = str_split(V5, "&", simplify = T)[,1],
         V8 = gsub("\\(.*\\)", "", V8), 
         V9 = gsub("\\(.*\\)", "", V9),
         PATHOGENIC = ifelse((V8 == "deleterious" | V9 == "probably_damaging") | V6 == "HIGH", "TRUE", NA)) %>%
  
  # filter variants to those that are high or moderate impact, are protein coding, and have a gene name
  # plus remove any variants in the mitochondrial genome 
  filter((V6 == "MODERATE" | V6 == "HIGH") & V7 == "protein_coding" & V4 != "NULL" & V4 != ".") %>%
  filter(!grepl("^MT-", V4)) %>%
  filter(V10 == "YES") %>%
  
  # filter to only distinct instances of each mutation
  dplyr::select(V11, V1, V2, V4, V5, PATHOGENIC) %>%
  distinct() %>%
  dplyr::select(V11, V4, V5, PATHOGENIC) %>%
  
  # rename columns and add in extra identifiers for oncoplot
  rename("SAMPLE" = "V11", "GENE" = "V4", "TYPE" = "V5") %>%
  mutate(CATEGORY = "Coding SNV/INDELs", CLUSTER = NA, ONCOPLOT = "SNV / SV")

##################
# Format SV data #
##################
svFormatted <- fread("data/svFormatted.bedpe") %>%
  dplyr::select(V1, V2, V3, V4, V5, V6, V14)

svClusters <- fread("data/Combined_consensus.100bp.clusters_foldbackinv_complex.txt") %>%
  dplyr::select(chr1, start1, end1, chr2, start2, end2, sv_type, clustered)

svAnnotated <- fread("data/svAnnotated.txt.gz")

svs <- svAnnotated %>%
  
  # annotate svs with breakpoints then annotate with clusters using breakpoints
  left_join(., svFormatted, by = c("id" = "V14")) %>%
  left_join(., svClusters, by = c("V1" = "chr1", "V2" = "start1", "V3" = "end1", "V4" = "chr2", "V5" = "start2", "V6" = "end2")) %>%
  mutate(CLUSTER = ifelse(clustered == 1, "Clustered SV", ifelse(clustered == 0, "Simple SV", NA)), 
         CATEGORY = "Structural Rearrangement", 
         PATHOGENIC = ifelse(exon_hit == "TRUE", "TRUE", NA)) %>%
  
  # filter to only svs disrupting exons and genes that are protein-coding
  filter(exon_hit == "TRUE" & biotype == "protein_coding") %>%
  dplyr::select(sample, gene, mut_type, PATHOGENIC, CATEGORY, CLUSTER) %>%
  
  # rename all columns and add in oncoplot identifier
  rename("SAMPLE" = "sample", "GENE" = "gene", "TYPE" = "mut_type") %>%
  mutate(ONCOPLOT = "SNV / SV")

###############################
# Format CNV duplication data #
###############################
cnvDupAnnotated <- fread("data/cnvDupAnnotated.txt.gz")

dups <- cnvDupAnnotated %>%

  # filter to only svs disrupting exons and genes that are protein-coding
  filter(exon_hit == "TRUE" & biotype == "protein_coding") %>%
  
  # add in extra identifiers and filter out any cnvs not hitting genes
  mutate(TYPE = "CNV_DUP", 
         CATEGORY = "CNV_DUP", 
         CLUSTER = NA, 
         PATHOGENIC = ifelse(exon_hit == "TRUE", "TRUE", NA)) %>%
  filter(gene != "NULL" & gene != ".") %>%
  
  # rename all columns and add in oncoplot identifier
  dplyr::select(sample, gene, TYPE, PATHOGENIC, CATEGORY, CLUSTER) %>%
  distinct() %>%
  rename("SAMPLE" = "sample", "GENE" = "gene") %>%
  mutate(ONCOPLOT = "CNV")

############################
# Format CNV deletion data #
############################
cnvDelAnnotated <- fread("data/cnvDelAnnotated.txt.gz")

dels <- cnvDelAnnotated %>%
  
  # filter to only svs disrupting exons and genes that are protein-coding
  filter(exon_hit == "TRUE" & biotype == "protein_coding") %>%
  
  # add in extra identifiers and filter out any cnvs not hitting genes
  mutate(TYPE = "CNV_DEL", CATEGORY = "CNV_DEL", 
         CLUSTER = NA, 
         PATHOGENIC = ifelse(exon_hit == "TRUE", "TRUE", NA)) %>%
  filter(gene != "NULL" & gene != ".") %>%
  
  # rename all columns and add in oncoplot identifier
  dplyr::select(sample, gene, TYPE, PATHOGENIC, CATEGORY, CLUSTER) %>%
  distinct() %>%
  rename("SAMPLE" = "sample", "GENE" = "gene") %>%
  mutate(ONCOPLOT = "CNV")

################
# Combine data #
################
combined <- do.call("rbind", list(snvs, svs, dups, dels))
dim(combined)
fwrite(combined, "data/oncoplotCombinedExonsHit.txt.gz", sep = "\t", quote = F, row.names = F)
