#####################################################
# Oncoplot grouped by pathways - Sanchez-Vega et al #
#####################################################

library(tidyverse)
library(ComplexHeatmap)
library(data.table)
library(GSA)


################
# SNV Pathways #
################
## Read in the Sanchez-Vega et al. gene set
gene_set <- GSA.read.gmt("data/sanchez_vega.gmt")

## Read in data, filter to pathogenic SNVs in genes in Sanchez-Vega et al. pathways 
## For each sample, make sure there is only one occurence of the gene, if there are multiple hits to a gene
## label these as "multi-hit"
snvs <- fread("data/oncoplotCombinedExonsHit.txt.gz") %>%
  filter(CATEGORY == "Coding SNV/INDELs" & PATHOGENIC == TRUE & GENE %in% unlist(gene_set$genesets)) %>%
  dplyr::select(SAMPLE, GENE, TYPE) %>%
  group_by(GENE, SAMPLE) %>%
  dplyr::mutate(TYPE = ifelse(length(SAMPLE) > 1, "multi_hit", TYPE)) %>%
  distinct()

## For each of the samples loop through eacg gene set and see if a gene within the set has been hit
pathways <- data.frame(pathway = gene_set$geneset.names)
for(i in unique(snvs$SAMPLE)){
  
  message(i)
  single <- snvs %>% filter(SAMPLE == i)
  
  sample <- data.frame()
  for(x in 1:length(gene_set$genesets)){
    
    genes <- gene_set$genesets[[x]]
    gene_set_name <- gene_set$geneset.names[x]
    
    ## use this to print what is hitting the genes in individual pathways
    single_pathway <- single %>% filter(GENE %in% genes)
    value <- ifelse(nrow(single_pathway) > 1, "multi_gene", single_pathway$TYPE)
    
    # value <- ifelse(sum(single$GENE %in% genes) != 0, paste0(gene_set_name,";"), "") # hit at least one gene
    # value <- ifelse(sum(single$GENE %in% genes) != 0, "HIT;", "") # hit at least one gene
    # value <- ifelse(sum(single$GENE %in% genes) > 1, "HIT;", "") # hit more than x genes
    # value <- ifelse((sum(single$GENE %in% genes) / length(genes)) * 100 > 5, "HIT;", "") # hit x% of genes
    
    out <- data.frame(pathway = gene_set_name, sample = value)
    colnames(out)[2] <- i
    
    sample <- rbind(sample, out)
    
  }
  
  pathways <- left_join(pathways, sample, by = "pathway")
  
}

# Format matrix of pathways
pathways <- pathways %>% column_to_rownames("pathway")

# Number of samples with at least one gene hit in each pathway
counts <- pathways %>%
  replace(. != "", 1) %>%
  mutate_all(as.numeric) %>%
  rowSums(na.rm = T) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  arrange(-.) %>%
  column_to_rownames("rowname") %>%
  head(n = 30)

pathways <- pathways %>% filter(row.names(pathways) %in% rownames(counts))

## Add in samples that do not have any pathway disrupted
usable <- fread("data/Combined_usable_ids_marker.txt", header = F) %>%
  mutate(V1 = gsub("_2$|-P$", "", V1))

pathways[setdiff(usable$V1, colnames(pathways))] <- ""
dim(pathways)

##### use this when using simple classification of pathways i.e. hit or not hit
#col = c("HIT" = "darkgreen")
# alter_fun <- list(
#   background = alter_graphic("rect", fill = "#F2F2F2", horiz_margin = unit(0.8, "points"), vertical_margin = unit(0.8, "points")),    
#   HIT = alter_graphic("rect", horiz_margin = unit(0.8, "points"), vertical_margin = unit(0.8, "points"), fill = col["HIT"])
# )

##### use this when using what is actually hitting individual pathways
col = c("frameshift_variant" = "#3BAE7C",
        "splice_donor_variant" = "#FF9FB2",
        "multi_hit" = "#242423",
        "missense_variant" = "#E25600",
        "stop_gained" = "#FCA910",
        "multi_gene" = "#6148B9",
        "splice_acceptor_variant" = "#737373", 
        "start_lost" = "#0000FF",
        "stop_lost" = "#FF6F0F")

## Format legend
heatmap_legend_param = list(title = "Alterations",
                            at = c("missense_variant", "frameshift_variant", "stop_gained",
                                   "multi_hit", "splice_acceptor_variant", "splice_donor_variant", "start_lost", "multi_gene"),
                            labels = c("Missense Variant", "Frameshift Variant", "Stop Gained",
                                       "Multi-Hit", "Splice-Acceptor Variant", "Splice-Donor Variant", "Start Lost", 'Multi-Gene-Hit'))

alter_fun <- list(
  background = alter_graphic("rect", fill = "#F2F2F2", horiz_margin = unit(0.3, "points"), vertical_margin = unit(0.3, "points")),
  frameshift_variant = alter_graphic("rect", horiz_margin = unit(0.3, "points"), vertical_margin = unit(0.3, "points"), fill = col["frameshift_variant"]),
  splice_donor_variant = alter_graphic("rect", horiz_margin = unit(0.3, "points"), vertical_margin = unit(0.3, "points"), fill = col["splice_donor_variant"]),
  multi_hit = alter_graphic("rect", horiz_margin = unit(0.3, "points"), vertical_margin = unit(0.3, "points"), fill = col["multi_hit"]),
  missense_variant = alter_graphic("rect", horiz_margin = unit(0.3, "points"), vertical_margin = unit(0.3, "points"), fill = col["missense_variant"]),
  stop_gained = alter_graphic("rect", horiz_margin = unit(0.3, "points"), vertical_margin = unit(0.3, "points"), fill = col["stop_gained"]),
  multi_gene = alter_graphic("rect", horiz_margin = unit(0.3, "points"), vertical_margin = unit(0.3, "points"), fill = col["multi_gene"]),
  splice_acceptor_variant = alter_graphic("rect", horiz_margin = unit(0.3, "points"), vertical_margin = unit(0.3, "points"), fill = col["splice_acceptor_variant"]),
  start_lost = alter_graphic("rect", horiz_margin = unit(0.3, "points"), vertical_margin = unit(0.3, "points"), fill = col["start_lost"]),
  stop_lost = alter_graphic("rect", horiz_margin = unit(0.3, "points"), vertical_margin = unit(0.3, "points"), fill = col["stop_lost"])
)

pdf(file = "snvsPathways.pdf", width = 20, height = 4.5)
oncoPrint(pathways,
          alter_fun = alter_fun, 
          col = col, 
          row_names_side = "left", 
          pct_side = "right", 
          pct_gp = gpar(fontsize = 16), 
          gap = unit(5, "mm"),
          heatmap_legend_param = heatmap_legend_param,
          top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(bar_width = 0.7),
                                             height = unit(2, units = "cm")),
          right_annotation = rowAnnotation(rbar = anno_oncoprint_barplot(bar_width = 0.8), 
                                           width = unit(3, units = "cm")))
dev.off()

###############
# SV Pathways #
###############
## Read in the Sanchez-Vega et al. gene set
gene_set <- GSA.read.gmt("data/sanchez_vega.gmt")

## Read in data, filter to pathogenic SNVs in genes in Sanchez-Vega et al. pathways 
## For each sample, make sure there is only one occurence of the gene, if there are multiple hits to a gene
## label these as "multi-hit"
svs <- fread("data/oncoplotCombinedExonsHit.txt.gz") %>%
  filter(CATEGORY == "Structural Rearrangement" & PATHOGENIC == TRUE & GENE %in% unlist(gene_set$genesets)) %>%
  dplyr::select(SAMPLE, GENE, TYPE) %>%
  group_by(GENE, SAMPLE) %>%
  dplyr::mutate(TYPE = ifelse(length(SAMPLE) > 1, "multi_hit", TYPE)) %>%
  distinct()

## For each of the samples loop through eacg gene set and see if a gene within the set has been hit
pathways <- data.frame(pathway = gene_set$geneset.names)
for(i in unique(svs$SAMPLE)){
  
  message(i)
  single <- svs %>% filter(SAMPLE == i)
  
  sample <- data.frame()
  for(x in 1:length(gene_set$genesets)){
    
    genes <- gene_set$genesets[[x]]
    gene_set_name <- gene_set$geneset.names[x]
    
    ## use this to print what is hitting the genes in individual pathways
    single_pathway <- single %>% filter(GENE %in% genes) 
    value <- ifelse(nrow(single_pathway) == 0, "", 
                    ifelse(nrow(single_pathway) == 1, single_pathway$TYPE,
                           ifelse(length(unique(single_pathway$TYPE)) == 1, paste0("multi_gene_", unique(single_pathway$TYPE)), "multi_gene")))
    
    # value <- ifelse(sum(single$GENE %in% genes) != 0, "HIT;", "") # hit at least one gene
    # value <- ifelse(sum(single$GENE %in% genes) >= 2, "HIT;", "") # hit more than x genes
    # value <- ifelse((sum(single$GENE %in% genes) / length(genes)) * 100 > 5, "HIT;", "") # hit x% of genes
    
    out <- data.frame(pathway = gene_set_name, sample = value)
    colnames(out)[2] <- i
    
    sample <- rbind(sample, out)
    
  }
  
  pathways <- left_join(pathways, sample, by = "pathway")
  
}

# Format matrix of pathways
pathways <- pathways %>% column_to_rownames("pathway")

## Add in samples that do not have any pathway disrupted
usable <- fread("data/Combined_usable_ids_marker.txt", header = F) %>%
  mutate(V1 = gsub("_2$|-P$", "", V1))

pathways[setdiff(usable$V1, colnames(pathways))] <- ""
dim(pathways)

# Number of samples with at least one gene hit in each pathway
counts <- pathways %>%
  replace(. != "", 1) %>%
  mutate_all(as.numeric) %>%
  rowSums(na.rm = T) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  arrange(-.) %>%
  column_to_rownames("rowname") %>%
  head(n = 30)

pathways <- pathways %>% filter(row.names(pathways) %in% rownames(counts))

#### use this when using simple classification of pathways i.e. hit or not hit
# col = c("HIT" = "darkgreen")
# alter_fun <- list(
#   background = alter_graphic("rect", fill = "#F2F2F2", horiz_margin = unit(0.8, "points"), vertical_margin = unit(0.8, "points")),
#   HIT = alter_graphic("rect", horiz_margin = unit(0.8, "points"), vertical_margin = unit(0.8, "points"), fill = col["HIT"])
# )

##### use this when using what is actually hitting individual pathways
col = c("multi_hit" = "#242423",
        "multi_gene" = "#6148B9",
        "DUP" = "#3AAF00",
        "INV" = "#A52A2A",
        "DEL" = "#00489F",
        "CTX" = "#D572BB",
        "multi_gene_DUP" = "#9BFF6A", 
        "multi_gene_DEL" = "#499EFF",
        "multi_gene_multi_hit" = "#A293D6", 
        "multi_gene_INV" = "#E08484")

## Format legend
heatmap_legend_param = list(title = "Alterations",
                            at = c("multi_hit", "multi_gene", "DUP", "INV", "DEL", "CTX", "multi_gene_DUP",  
                                   "multi_gene_DEL", "multi_gene_multi_hit", "multi_gene_INV"),
                            labels = c("Multi-Hit", "Multi-Gene-Hit", "Duplication", "Inversion", "Deletion", "Translocation", "Multi-Gene-Duplication",  
                                       "Multi-Gene-Deletion", "Multi-Gene-Multi-Hit", "Multi-Gene-Inversion"))

alter_fun <- list(
  background = alter_graphic("rect", fill = "#F2F2F2", horiz_margin = unit(0.3, "points"), vertical_margin = unit(0.3, "points")),
  multi_hit = alter_graphic("rect", horiz_margin = unit(0.3, "points"), vertical_margin = unit(0.3, "points"), fill = col["multi_hit"]),
  multi_gene = alter_graphic("rect", horiz_margin = unit(0.3, "points"), vertical_margin = unit(0.3, "points"), fill = col["multi_gene"]),
  DUP = alter_graphic("rect", horiz_margin = unit(0.3, "points"), vertical_margin = unit(0.3, "points"), fill = col["DUP"]),
  INV = alter_graphic("rect", horiz_margin = unit(0.3, "points"), vertical_margin = unit(0.3, "points"), fill = col["INV"]),
  DEL = alter_graphic("rect", horiz_margin = unit(0.3, "points"), vertical_margin = unit(0.3, "points"), fill = col["DEL"]),
  CTX = alter_graphic("rect", horiz_margin = unit(0.3, "points"), vertical_margin = unit(0.3, "points"), fill = col["CTX"]),
  multi_gene_DUP = alter_graphic("rect", horiz_margin = unit(0.3, "points"), vertical_margin = unit(0.3, "points"), fill = col["multi_gene_DUP"]),
  multi_gene_DEL = alter_graphic("rect", horiz_margin = unit(0.3, "points"), vertical_margin = unit(0.3, "points"), fill = col["multi_gene_DEL"]),
  multi_gene_multi_hit = alter_graphic("rect", horiz_margin = unit(0.3, "points"), vertical_margin = unit(0.3, "points"), fill = col["multi_gene_multi_hit"]),
  multi_gene_INV = alter_graphic("rect", horiz_margin = unit(0.3, "points"), vertical_margin = unit(0.3, "points"), fill = col["multi_gene_INV"])
)

pdf(file = "svsPathways.pdf", width = 20, height = 4.5)
oncoPrint(pathways,
          alter_fun = alter_fun, 
          col = col, 
          row_names_side = "left", 
          pct_side = "right", 
          pct_gp = gpar(fontsize = 16), 
          gap = unit(5, "mm"),
          heatmap_legend_param = heatmap_legend_param,
          top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(bar_width = 0.7),
                                             height = unit(2, units = "cm")),
          right_annotation = rowAnnotation(rbar = anno_oncoprint_barplot(bar_width = 0.8), 
                                           width = unit(3, units = "cm")))
dev.off()


################
# CNV Pathways #
################
## Read in the Sanchez-Vega et al. gene set
gene_set <- GSA.read.gmt("data/sanchez_vega.gmt")

## Read in data, filter to pathogenic SNVs in genes in Sanchez-Vega et al. pathways 
## For each sample, make sure there is only one occurence of the gene, if there are multiple hits to a gene
## label these as "multi-hit"
cnvs <- fread("data/oncoplotCombinedExonsHit.txt.gz") %>%
  filter(CATEGORY %in% c("CNV_DUP", "CNV_DEL") & PATHOGENIC == TRUE & GENE %in% unlist(gene_set$genesets)) %>%
  dplyr::select(SAMPLE, GENE, TYPE) %>%
  group_by(GENE, SAMPLE) %>%
  dplyr::mutate(TYPE = ifelse(length(SAMPLE) > 1, "multi_hit", TYPE)) %>%
  distinct()

## For each of the samples loop through eacg gene set and see if a gene within the set has been hit
pathways <- data.frame(pathway = gene_set$geneset.names)
for(i in unique(cnvs$SAMPLE)){
  
  message(i)
  single <- cnvs %>% filter(SAMPLE == i)
  
  sample <- data.frame()
  for(x in 1:length(gene_set$genesets)){
    
    genes <- gene_set$genesets[[x]]
    gene_set_name <- gene_set$geneset.names[x]
    
    ## use this to print what is hitting the genes in individual pathways
    single_pathway <- single %>% filter(GENE %in% genes) 
    value <- ifelse(nrow(single_pathway) == 0, "", 
                    ifelse(nrow(single_pathway) == 1, single_pathway$TYPE,
                           ifelse(length(unique(single_pathway$TYPE)) == 1, paste0("multi_gene_", unique(single_pathway$TYPE)), "multi_gene")))
    
    # value <- ifelse(sum(single$GENE %in% genes) != 0, "HIT;", "") # hit at least one gene
    # value <- ifelse(sum(single$GENE %in% genes) >= 5, "HIT;", "") # hit more than x genes
    # value <- ifelse((sum(single$GENE %in% genes) / length(genes)) * 100 > 10, "HIT;", "") # hit x% of genes
    
    out <- data.frame(pathway = gene_set_name, sample = value)
    colnames(out)[2] <- i
    
    sample <- rbind(sample, out)
    
  }
  
  pathways <- left_join(pathways, sample, by = "pathway")
  
}

# Format matrix of pathways
pathways <- pathways %>% column_to_rownames("pathway")

## Add in samples that do not have any pathway disrupted
usable <- fread("data/Combined_usable_ids_marker.txt", header = F) %>%
  mutate(V1 = gsub("_2$|-P$", "", V1))

pathways[setdiff(usable$V1, colnames(pathways))] <- ""
dim(pathways)

# Number of samples with at least one gene hit in each pathway
counts <- pathways %>%
  replace(. != "", 1) %>%
  mutate_all(as.numeric) %>%
  rowSums(na.rm = T) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  arrange(-.) %>%
  column_to_rownames("rowname") %>%
  head(n = 20)

pathways <- pathways %>% filter(row.names(pathways) %in% rownames(counts))

#### use this when using simple classification of pathways i.e. hit or not hit
# col = c("HIT" = "darkgreen")
# alter_fun <- list(
#   background = alter_graphic("rect", fill = "#F2F2F2", horiz_margin = unit(0.8, "points"), vertical_margin = unit(0.8, "points")),
#   HIT = alter_graphic("rect", horiz_margin = unit(0.8, "points"), vertical_margin = unit(0.8, "points"), fill = col["HIT"])
# )

col = c("multi_gene_CNV_DUP" = "#9BFF6A", 
        "CNV_DUP" = "#3AAF00",
        "multi_gene_CNV_DEL" = "#499EFF", 
        "multi_gene" = "#6148B9", 
        "CNV_DEL" = "#00489F")

## Format legend
heatmap_legend_param = list(title = "Alterations",
                            at = c("multi_gene_CNV_DUP", "CNV_DUP", "multi_gene_CNV_DEL", "multi_gene", "CNV_DEL"),
                            labels = c("Multi-Gene-Duplication", "Duplication", "Multi-Gene-Deletion", "Multi-Gene", "Deletion"))

alter_fun <- list(
  background = alter_graphic("rect", fill = "#F2F2F2", horiz_margin = unit(0.3, "points"), vertical_margin = unit(0.3, "points")),
  multi_gene_CNV_DUP = alter_graphic("rect", horiz_margin = unit(0.3, "points"), vertical_margin = unit(0.3, "points"), fill = col["multi_gene_CNV_DUP"]),
  CNV_DUP = alter_graphic("rect", horiz_margin = unit(0.3, "points"), vertical_margin = unit(0.3, "points"), fill = col["CNV_DUP"]),
  multi_gene_CNV_DEL = alter_graphic("rect", horiz_margin = unit(0.3, "points"), vertical_margin = unit(0.3, "points"), fill = col["multi_gene_CNV_DEL"]),
  multi_gene = alter_graphic("rect", horiz_margin = unit(0.3, "points"), vertical_margin = unit(0.3, "points"), fill = col["multi_gene"]),
  CNV_DEL = alter_graphic("rect", horiz_margin = unit(0.3, "points"), vertical_margin = unit(0.3, "points"), fill = col["CNV_DEL"])
)

pdf(file = "cnvsPathways.pdf", width = 20, height = 4.5)
oncoPrint(pathways,
          alter_fun = alter_fun, 
          col = col, 
          row_names_side = "left", 
          pct_side = "right", 
          pct_gp = gpar(fontsize = 16), 
          gap = unit(5, "mm"),
          heatmap_legend_param = heatmap_legend_param,
          top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(bar_width = 0.7),
                                             height = unit(2, units = "cm")),
          right_annotation = rowAnnotation(rbar = anno_oncoprint_barplot(bar_width = 0.8), 
                                           width = unit(3, units = "cm")))
dev.off()
