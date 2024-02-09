## Figure 1##

##########################
# Introduce cohort plots #
##########################

library(data.table)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(gtools)
library(purrr)
library(scales)
library(tibble)
library(tidyr)



##### Usable Samples #####
usable <- fread("data/Combined_usable_ids_marker.txt", header = F) %>%
  mutate(V1 = gsub("_2$|-P$", "", V1))

##### Format cSV Matrix #####
csv <- fread("data/SummaryTableSampleCallOfcSV.tsv", header = T) %>%
  dplyr::select(-V1, -Cohort, -dplyr::contains("Starfish")) %>%
  mutate(Sample = gsub("-P$", "", Sample), SS_Chromothripsis = ifelse(SS_Chromothripsis != "Absent", "Present", SS_Chromothripsis)) 

##### Format Expression Subtypes Matrix #####
exp <- fread("data/expression_subtypes.csv") %>%
  dplyr::select(sample_id, decision.consensus) %>%
  mutate(sample_id = gsub("T$|-P$", "", sample_id)) %>%
  add_row(sample_id = setdiff(usable$V1, .[["sample_id"]]), decision.consensus = "UNKNOWN") %>%
  dplyr::rename("Sample" = "sample_id")

##### Format Cohort and Stage Matrix #####
coh <- fread("data/HGSOC_cohorts_minimal_clin.txt") %>%
  dplyr::select(Sample, stage, cohort) %>%
  mutate(cohort = ifelse(cohort == "SHGSOC_new", "SHGSOC", cohort),
         Sample = gsub("-P$|_2$", "", Sample), stage = paste0("stage", stage))

##### Format Germline SNVs
germline <- fread("data/germlineKittyMutations.txt") %>%
  mutate(Sample = gsub("N$|-[0-9]$", "", Sample)) %>%
  dplyr::select(Sample, Gene) %>% 
  distinct() %>%
  aggregate(data = ., .~Sample, FUN = paste, collapse = ";")

##### Combine Matrices 
list_complex <- list(csv, exp, coh, germline)
combined <- list_complex %>% 
  purrr::reduce(full_join, by = "Sample") %>%
  column_to_rownames("Sample") %>% t() %>%
  replace(is.na(.), "")

order <- c("AA_ecDNA" = "A", "SS_Chromothripsis" = "A", "WholeGenomeDuplication" = "B",
           "rigma" = "A", "pyrgo" = "A", "chromoplexy" = "A", "bfb" = "A", "tyfonas" = "A",
           "HRDetectStat" = "B", "SeismicAmplification" = "A", 
           #"decision.consensus" = "C",
           "stage" = "B", "cohort" = "B", "Gene" = "B")

##### Tidy up Row Names #####
rownames(combined)[rownames(combined) == "AA_ecDNA"] <- "ecDNA"
rownames(combined)[rownames(combined) == "SS_Chromothripsis"] <- "Chromothripsis"
rownames(combined)[rownames(combined) == "WholeGenomeDuplication"] <- "WGD"
rownames(combined)[rownames(combined) == "rigma"] <- "Rigma"
rownames(combined)[rownames(combined) == "pyrgo"] <- "Pyrgo"
rownames(combined)[rownames(combined) == "chromoplexy"] <- "Chromoplexy"
rownames(combined)[rownames(combined) == "bfb"] <- "BFB"
rownames(combined)[rownames(combined) == "tyfonas"] <- "Tyfonas"
rownames(combined)[rownames(combined) == "HRDetectStat"] <- "HRD Status"
rownames(combined)[rownames(combined) == "SeismicAmplification"] <- "Seismic Amplification"
rownames(combined)[rownames(combined) == "decision.consensus"] <- "Expression Group"
rownames(combined)[rownames(combined) == "stage"] <- "Stage"
rownames(combined)[rownames(combined) == "cohort"] <- "Cohort"
rownames(combined)[rownames(combined) == "Gene"] <- "Germline"

##### Define Plotting Information for Oncoplot #####


col = c("Present" = "#242423",
        "Absent" = "#F2F2F2",
      #  "MES" = "#9b5de5", 
      #  "DIF" = "#f15bb5", 
      #  "PRO" = "#fee440",
      #  "IMR" = "#00bbf9",
      #  "UNKNOWN" = "#F2F2F2",
        "AOCS" = "blue",
        "BCCA" = "black",
        "SHGSOC" = "brown",
        "TCGA" = "grey",
        "MDA" = "orange",
        "stage1" = "#FEE5D9",
        "stage2" = "#FCAE91", 
        "stage3" = "#FB6A4A", 
        "stage4" = "#7F000D", 
        "stageNA" = "#F2F2F2",
        "BRCA1" = "#F44336",
        "BRCA2" = "#242423",
        "BRIP1" = "#1565c0",
        "ATM" = "#009688",
        "MSH6" = "#8bc34a",
        "RAD51C" = "#ffc107",
        "RAD51D" = "#ff9800",
        "PMS2" = "#ad1457", 
        "PALB2" = "dodgerblue")

alter_fun <- list(
  
  background = alter_graphic("rect", fill = "#F2F2F2", horiz_margin = unit(0.2, "points"), vertical_margin = unit(0.2, "points")),   
  Absent = alter_graphic("rect", horiz_margin = unit(0.2, "points"), vertical_margin = unit(0.2, "points"), fill = col["Absent"]),
  Present = alter_graphic("rect", horiz_margin = unit(0.2, "points"), vertical_margin = unit(0.2, "points"), fill = col["Present"]),
  
  #MES = alter_graphic("rect", horiz_margin = unit(0.2, "points"), vertical_margin = unit(0.2, "points"), fill = col["MES"]),
 # DIF = alter_graphic("rect", horiz_margin = unit(0.2, "points"), vertical_margin = unit(0.2, "points"), fill = col["DIF"]),
  #PRO = alter_graphic("rect", horiz_margin = unit(0.2, "points"), vertical_margin = unit(0.2, "points"), fill = col["PRO"]),
  #IMR = alter_graphic("rect", horiz_margin = unit(0.2, "points"), vertical_margin = unit(0.2, "points"), fill = col["IMR"]),
  #UNKNOWN = alter_graphic("rect", horiz_margin = unit(0.2, "points"), vertical_margin = unit(0.2, "points"), fill = col["UNKNOWN"]),
  
  AOCS = alter_graphic("rect", horiz_margin = unit(0.2, "points"), vertical_margin = unit(0.2, "points"), fill = col["AOCS"]),
  BCCA = alter_graphic("rect", horiz_margin = unit(0.2, "points"), vertical_margin = unit(0.2, "points"), fill = col["BCCA"]),
  SHGSOC = alter_graphic("rect", horiz_margin = unit(0.2, "points"), vertical_margin = unit(0.2, "points"), fill = col["SHGSOC"]),
  TCGA = alter_graphic("rect", horiz_margin = unit(0.2, "points"), vertical_margin = unit(0.2, "points"), fill = col["TCGA"]),
  MDA = alter_graphic("rect", horiz_margin = unit(0.2, "points"), vertical_margin = unit(0.2, "points"), fill = col["MDA"]),
  
  stage1 = alter_graphic("rect", horiz_margin = unit(0.2, "points"), vertical_margin = unit(0.2, "points"), fill = col["stage1"]),
  stage2 = alter_graphic("rect", horiz_margin = unit(0.2, "points"), vertical_margin = unit(0.2, "points"), fill = col["stage2"]),
  stage3 = alter_graphic("rect", horiz_margin = unit(0.2, "points"), vertical_margin = unit(0.2, "points"), fill = col["stage3"]),
  stage4 = alter_graphic("rect", horiz_margin = unit(0.2, "points"), vertical_margin = unit(0.2, "points"), fill = col["stage4"]),
  stageNA = alter_graphic("rect", horiz_margin = unit(0.2, "points"), vertical_margin = unit(0.2, "points"), fill = col["stageNA"]),
  
  BRCA1 = alter_graphic("rect", horiz_margin = unit(0.2, "points"), vertical_margin = unit(0.2, "points"), fill = col["BRCA1"]),
  BRIP1 = alter_graphic("rect", horiz_margin = unit(0.2, "points"), vertical_margin = unit(0.2, "points"), fill = col["BRIP1"]),
  BRCA2 = alter_graphic("rect", horiz_margin = unit(0.2, "points"), vertical_margin = unit(0.2, "points"), fill = col["BRCA2"]),
  ATM = alter_graphic("rect", horiz_margin = unit(0.2, "points"), vertical_margin = unit(0.2, "points"), fill = col["ATM"]),
  MSH6 = alter_graphic("rect", horiz_margin = unit(0.2, "points"), vertical_margin = unit(0.2, "points"), fill = col["MSH6"]),
  RAD51C = alter_graphic("rect", horiz_margin = unit(0.2, "points"), vertical_margin = unit(0.2, "points"), fill = col["RAD51C"]),
  RAD51D = alter_graphic("rect", horiz_margin = unit(0.2, "points"), vertical_margin = unit(0.2, "points"), fill = col["RAD51D"]),
  PMS2 = alter_graphic("rect", horiz_margin = unit(0.2, "points"), vertical_margin = unit(0.2, "points"), fill = col["PMS2"]),
  PALB2 = alter_graphic("rect", horiz_margin = unit(0.2, "points"), vertical_margin = unit(0.2, "points"), fill = col["PALB2"])
  
)

##### Create Top Annotation #####
snvs <- fread("data/simpleSNVs.txt.gz") %>%
  mutate(V4 = gsub("-ensemble-annotated-oxog_filtered.vep107.vcf.gz", "", V4)) %>%
  mutate(V4 = gsub("_2$|-P$", "", V4)) %>%
  filter(V1 != "chrM") %>%
  group_by(V4) %>%
  tally() %>%
  dplyr::rename("sample" = "V4", "SNVs" = "n")

svs <- fread("data/svFormatted.bedpe") %>%
  mutate(V12 = abs(V12)) %>%
  mutate(V13 = gsub("_PrimaryTumour", "", V13)) %>%
  mutate(V13 = gsub("T$", "", V13)) %>% 
  mutate(V13 = gsub("T2$", "", V13)) %>%
  mutate(V13 = gsub("-[0-9]$", "", V13)) %>%
  dplyr::select(V13, V12, V11) %>%
  aggregate(V12 ~ V13 + V11, data = ., sum) %>%
  pivot_wider(names_from = V11, values_from = V12) %>%
  dplyr::rename("sample" = "V13")

translocations <- fread("data/svFormatted.bedpe") %>%
  mutate(V13 = gsub("_PrimaryTumour", "", V13)) %>%
  mutate(V13 = gsub("T$", "", V13)) %>% 
  mutate(V13 = gsub("T2$", "", V13)) %>%
  mutate(V13 = gsub("-[0-9]$", "", V13)) %>%
  dplyr::select(V13, V12, V11) %>%
  filter(V11 == "CTX") %>%
  group_by(V13, V11) %>%
  tally() %>%
  pivot_wider(names_from = V11, values_from = n) %>%
  dplyr::rename("sample" = "V13")

dels <- fread("data/delFormatted.bed") %>%
  mutate(bases = V3 - V2,
         mut_type = "CNV_DEL") %>%
  dplyr::select(V5, bases, mut_type) %>%
  aggregate(bases ~ V5 + mut_type, data = ., sum) %>%
  pivot_wider(names_from = mut_type, values_from = bases) %>%
  dplyr::rename("sample" = "V5")

dups <- fread("data/dupFormatted.bed") %>%
  mutate(bases = V3 - V2,
         mut_type = "CNV_DUP") %>%
  dplyr::select(V5, bases, mut_type) %>%
  aggregate(bases ~ V5 + mut_type, data = ., sum) %>%
  pivot_wider(names_from = mut_type, values_from = bases) %>%
  dplyr::rename("sample" = "V5")


list_muts <- list(snvs, svs, translocations, dels, dups)
mutation_counts <- list_muts %>% purrr::reduce(full_join, by = "sample")
mutation_counts[is.na(mutation_counts)] <- 0

other_counts <- mutation_counts %>% dplyr::select(-SNVs, -CTX) %>% column_to_rownames("sample")
snv_counts <- mutation_counts %>% dplyr::select(sample, CTX) %>% column_to_rownames("sample")

ordered <- other_counts %>%
  mutate(total = rowSums(.)) %>%
  rownames_to_column("rownames") %>%
  arrange(total) %>%
  column_to_rownames("rownames") %>%
  dplyr::select(-total)

combined <- combined[,match(rownames(ordered), colnames(combined))]
snv_counts <- snv_counts[match(rownames(ordered), rownames(snv_counts)),]

barplot <- HeatmapAnnotation(other = anno_barplot(ordered, 
                                                  gp = gpar(fill = c(	"#D8BFD8", "#CC79A7", "#3AAF00","#ffc107", "#0072B2","#009E73"),
                                                          col = NA),
                                                  bar_width = 0.8, 
                                                  #ylim = c(0, 2500), 
                                                  border = F,
                                                  height = unit((11/3) * 2, "cm"),
                                                  axis_param = list(at=c(0,2e9,4e9,6e9,8e9),
                                                                    labels = c("0","2,000,000,000","4,000,000,000","6,000,000,000","8,000,000,000")
                                                                    ),
                                                  cex.axis=2
                                                  ),
                             snvs = anno_barplot(snv_counts, 
                                                 gp = gpar(fill = c("#242423"), 
                                                           col = NA), 
                                                 bar_width = 0.8, 
                                                 #ylim = c(0, 2500), 
                                                 #axis_param = list(at = seq(0, 2500, 500)), 
                                                 border = F,
                                                 height = unit(11/3, "cm")),
                             height = unit(11, units = "cm"),
                             show_legend = c(F,F),
                             show_annotation_name=FALSE)

###### Create Annotation For Deconvoluted Cellular Proportions ######
cells <- fread("data/Deconvolution_results.txt") %>%
  dplyr::select(-P.value, -Correlation, -RMSE) %>%
  mutate(Mixture = gsub("-P$", "", Mixture)) %>%
  column_to_rownames("Mixture") 

#cells[] <- t(apply(cells, 1, function(x){ x / sum(x) }))
cells <- cells[match(rownames(ordered), rownames(cells)),]
cells$unknown <- ifelse(is.na(cells$Endothelial), 1, 0)
cells[is.na(cells)] <- 0
cells <- cells %>% dplyr::rename("VSMC" = "VSMCs")

cell_proportions <- HeatmapAnnotation(foo = anno_barplot(cells, 
                                                         gp = gpar(fill = c("#F94144", "#1E90FF", "#F8961E", "#AD1457", "#F9C74F", "#90BE6D", "#43AA8B", "#F2F2F2"), 
                                                                   col = NA), 
                                                         bar_width = 0.8, 
                                                         border = F), 
                                      height = unit(2, units = "cm"),
                                      show_legend = T)


row_order <- c("ecDNA", "Chromothripsis", "Rigma", "Pyrgo", 'Chromoplexy', "BFB", "Tyfonas", "Seismic Amplification",
               "WGD", "HRD Status", "Germline", "Stage", "Cohort")
              # "Expression Group")


##### Format Legend #####
lgd1_cols = c("Deletion" = "#D8BFD8", "Duplication" = "#CC79A7", "Insertion" = "#3AAF00", "Inversion" = "#ffc107", "CNA Deletion" = "#0072B2", "CNA Duplication" = "#009E73")
lgd2_cols = c("Present" = "#242423", "Absent" = "#F2F2F2")
lgd3_cols = c("BRCA1" = "#F44336", "BRCA2" = "#242423", "BRIP1" = "#1565c0", "ATM" = "#009688", "MSH6" = "#8bc34a", "RAD51C" = "#ffc107", "RAD51D" = "#ff9800", "PMS2" = "#ad1457", "PALB2" = "dodgerblue")
lgd4_cols = c("Stage I" = "#FEE5D9","Stage II" = "#FCAE91", "Stage III" = "#FB6A4A", "Stage IV" = "#7F000D", "Unknown" = "#F2F2F2")
lgd5_cols = c("AOCS" = "blue", "BCCA" = "black", "SHGSOC" = "brown", "TCGA" = "grey", "MDA" = "orange")
#lgd6_cols = c("MES" = "#9b5de5", "DIF" = "#f15bb5", "PRO" = "#fee440", "IMR" = "#00bbf9", "No RNA-Seq" = "#F2F2F2")
#lgd7_cols = c("Endothelial" = "#F94144", "Epithelial" = "#1E90FF", "Fibroblast" = "#F8961E", "Stromal" = "#AD1457", "VSMC" = "#F9C74F", "Myeloid" = "#90BE6D", "Tcell" = "#43AA8B", "Unknown" = "#F2F2F2")

lgd1 = Legend(at = names(lgd1_cols), legend_gp = gpar(fill = lgd1_cols), labels_gp = gpar(fontsize=12),title = "SVs / CNAs",title_gp = gpar(fontsize = 12, fontface = "bold"))
lgd2 = Legend(at = names(lgd2_cols), legend_gp = gpar(fill = lgd2_cols), labels_gp = gpar(fontsize=12),title = "Complex Events",title_gp = gpar(fontsize = 12, fontface = "bold"))
lgd3 = Legend(at = names(lgd3_cols), legend_gp = gpar(fill = lgd3_cols), labels_gp = gpar(fontsize=12), title = "Germline",title_gp = gpar(fontsize = 12, fontface = "bold"))
lgd4 = Legend(at = names(lgd4_cols), legend_gp = gpar(fill = lgd4_cols), labels_gp = gpar(fontsize=12),title = "Tumour Stage",title_gp = gpar(fontsize = 12, fontface = "bold"))
lgd5 = Legend(at = names(lgd5_cols), legend_gp = gpar(fill = lgd5_cols), labels_gp = gpar(fontsize=12),title = "Cohort",title_gp = gpar(fontsize = 12, fontface = "bold"))
#lgd6 = Legend(at = names(lgd6_cols), legend_gp = gpar(fill = lgd6_cols), title = "Expression Group")
#lgd7 = Legend(at = names(lgd7_cols), legend_gp = gpar(fill = lgd7_cols), title = "Cell Types")

pd = packLegend(list = list(lgd1, lgd2, lgd3, lgd4, lgd5), direction = "vertical")

##### Plot Oncoplot #####

oncoplot = oncoPrint(combined[setdiff(rownames(combined),"Expression Group"),],
          alter_fun = alter_fun, 
          col = col, 
          row_names_side = "left", 
          pct_side = "right", 
          pct_gp = gpar(fontsize = 13), 
          show_pct = F,
          gap = unit(5, "mm"),
          row_split = order,
          row_title = NULL,
          top_annotation = barplot, 
      #   bottom_annotation = cell_proportions,
          column_order = rownames(ordered),
          right_annotation = NULL, 
          row_order = row_order)

pdf(file = "individual_panels/main_figs/Fig1/Fig1A.pdf", width = 20, height = 7)
draw(oncoplot, ht_gap = unit(7, "mm"), row_km = 2, show_heatmap_legend = FALSE, annotation_legend_list = pd)
dev.off()



