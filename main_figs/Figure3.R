###################################
# Driver Figure - Panel A / B / C #
###################################

library(ComplexHeatmap)
library(cowplot)
library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

#############################
# Driver Oncoplot - Panel A #
#############################
muts <- fread("data/oncoplotCombinedExonsHit.txt.gz")

sv_drivers <- c("CEP89", "CCNE1", "PTEN")
snv_drivers <- c("TP53", "NF1", "RB1", "SLC35G5", "BRCA1", "CDK12", "TAS2R43", "BRCA2") 
cnv_drivers <- c("AXIN1", "PCM1", "WRN", "CEP43", "LEPROTL1", "DNMT3A", "ARID1B", "TCEA1")

drivers <- unique(c(sv_drivers, snv_drivers, cnv_drivers))

driver_muts <- filter(muts, ((CATEGORY == "Structural Rearrangement" | CATEGORY == "Coding SNV/INDELs") & GENE %in% sv_drivers) |
                        ((CATEGORY == "Structural Rearrangement" | CATEGORY == "Coding SNV/INDELs") & GENE %in% snv_drivers) |
                        ((CATEGORY == "CNV_DEL" | CATEGORY == "CNV_DUP") & GENE %in% cnv_drivers))

matrix <- driver_muts %>%
  group_by(SAMPLE, GENE, CATEGORY) %>%
  dplyr::mutate(TYPE = ifelse(length(SAMPLE) > 1, "multi_hit", TYPE),
                TYPE = ifelse(TYPE == "multi_hit" & CATEGORY == "Structural Rearrangement", "multi_hit_structural_variation", TYPE)) %>%
  ungroup() %>% dplyr::select(-CATEGORY, -CLUSTER, -ONCOPLOT, -PATHOGENIC) %>%
  distinct() %>%
  aggregate(data = ., .~SAMPLE+GENE, FUN = paste, collapse = ";") %>%
  pivot_wider(names_from = SAMPLE, values_from = TYPE) %>%
  column_to_rownames("GENE")

matrix[] <- lapply(matrix, function(x){ ifelse(is.na(x), "", paste0(x, ";")) })

usable <- fread("data/Combined_usable_ids_marker.txt", header = F) %>%
  mutate(V1 = gsub("_2$|-P$", "", V1))

matrix[setdiff(usable$V1, colnames(matrix))] <- ""
dim(matrix)

subplots <- data.frame(gene = rownames(matrix)) %>%
  mutate(sub_plot = ifelse(gene %in% cnv_drivers, "CNV", "SNV / SV"),
         sub_plot = factor(sub_plot, levels=c("SNV / SV", "CNV"))) %>%
  arrange(sub_plot)

matrix <- matrix[match(subplots$gene, rownames(matrix)),]

##### Calculate Percentages #####
deletions <- c("TP53", "NF1", "PTEN", "RB1", "BRCA1", "CDK12", "TAS2R43", "SLC35G5", "BRCA2")
duplications <- c("CCNE1", "CEP89")
gain <- c("DNMT3A", "TCEA1")
loss <- c("LEPROTL1", "CEP43", "AXIN1", "WRN", "PCM1", "ARID1B")

pathogenic <- driver_muts %>%
  filter((GENE %in% deletions & ((CATEGORY == "Coding SNV/INDELs" & PATHOGENIC == "TRUE") | TYPE == "DEL")) |
         (GENE %in% duplications & ((CATEGORY == "Coding SNV/INDELs" & PATHOGENIC == "TRUE") | TYPE == "DUP")) |
         (GENE %in% gain & TYPE == "CNV_DUP") | 
         (GENE %in% loss & TYPE == "CNV_DEL")) %>%
  dplyr::select(SAMPLE, GENE) %>% distinct() %>%
  group_by(GENE) %>%
  tally() %>%
  mutate(percentage = round((n / 324) * 100))

overall <- driver_muts %>%
  dplyr::select(SAMPLE, GENE) %>% distinct() %>%
  group_by(GENE) %>%
  tally() %>%
  mutate(percentage = round((n / 324) * 100))

pathogenic <- pathogenic[match(rownames(matrix), pathogenic$GENE),]
overall <- overall[match(rownames(matrix), overall$GENE),]

col = c("missense_variant" = "#137547", 
        "frameshift_variant" = "#0B6FE3", 
        "stop_gained" = "#FBAC11", 
        "multi_hit" = "#242423", 
        "splice_acceptor_variant" = "#FF021B", 
        "inframe_deletion" = "#FFAB86",
        "inframe_insertion" = "#380099", 
        "splice_donor_variant" = "#84714F",
        
        "DUP" = "#CC79A7", ##
        "DEL" = "#D8BFD8", ##
        "INV" = "#56B4E9", 
        "CTX" = "#E69F00",
        
        "CNV_DUP" = "#009E73", ##
        "CNV_DEL" = "#0072B2", ##
        "multi_hit_structural_variation" = "#242423") ##


alter_fun <- list(
  
  ## Background colour
  background = alter_graphic("rect", fill = "#F2F2F2", horiz_margin = unit(0.3, "points"), vertical_margin = unit(0.3, "points")),   
  
  ## Somatic SNVs 
  missense_variant = alter_graphic("rect", horiz_margin = unit(0.3, "points"), vertical_margin = unit(0.3, "points"), fill = col["missense_variant"]),
  frameshift_variant = alter_graphic("rect", horiz_margin = unit(0.3, "points"), vertical_margin = unit(0.3, "points"), fill = col["frameshift_variant"]),
  stop_gained = alter_graphic("rect", horiz_margin = unit(0.3, "points"), vertical_margin = unit(0.3, "points"), fill = col["stop_gained"]),
  multi_hit = alter_graphic("rect", horiz_margin = unit(0.3, "points"), vertical_margin = unit(0.3, "points"), fill = col["multi_hit"]),
  splice_acceptor_variant = alter_graphic("rect", horiz_margin = unit(0.3, "points"), vertical_margin = unit(0.3, "points"), fill = col["splice_acceptor_variant"]),
  inframe_deletion = alter_graphic("rect", horiz_margin = unit(0.3, "points"), vertical_margin = unit(0.3, "points"), fill = col["inframe_deletion"]),
  inframe_insertion = alter_graphic("rect", horiz_margin = unit(0.3, "points"), vertical_margin = unit(0.3, "points"), fill = col["inframe_insertion"]),
  splice_donor_variant = alter_graphic("rect", horiz_margin = unit(0.3, "points"), vertical_margin = unit(0.3, "points"), fill = col["splice_donor_variant"]),
  
  ## SVs
  DUP = alter_graphic("rect", horiz_margin = unit(0.3, "points"), vertical_margin = unit(0.3, "points"), height = 0.5, fill = col["DUP"]),
  CTX = alter_graphic("rect", horiz_margin = unit(0.3, "points"), vertical_margin = unit(0.3, "points"), height = 0.5, fill = col["CTX"]),
  DEL = alter_graphic("rect", horiz_margin = unit(0.3, "points"), vertical_margin = unit(0.3, "points"), height = 0.5, fill = col["DEL"]),
  INV = alter_graphic("rect", horiz_margin = unit(0.3, "points"), vertical_margin = unit(0.3, "points"), height = 0.5, fill = col["INV"]),
  INS = alter_graphic("rect", horiz_margin = unit(0.3, "points"), vertical_margin = unit(0.3, "points"), height = 0.5, fill = col["INS"]),
  BND = alter_graphic("rect", horiz_margin = unit(0.3, "points"), vertical_margin = unit(0.3, "points"), height = 0.5, fill = col["BND"]),
  multi_hit_structural_variation = alter_graphic("rect", horiz_margin = unit(0.3, "points"), vertical_margin = unit(0.3, "points"), height = 0.5, fill = col["multi_hit"]),
  
  ## CNVs 
  CNV_DUP = alter_graphic("rect", horiz_margin = unit(0.3, "points"), vertical_margin = unit(0.3, "points"), height = 0.5, fill = col["CNV_DUP"]),
  CNV_DEL = alter_graphic("rect", horiz_margin = unit(0.3, "points"), vertical_margin = unit(0.3, "points"), height = 0.5, fill = col["CNV_DEL"])
  
)

##### Format Legend #####
lgd1_cols <- c("Missense" = "#137547", "Frameshift" = "#0B6FE3", "Stop Gained" = "#FBAC11", "Splice-Acceptor" = "#FF021B", "Inframe Deletion" = "#FFAB86",
          "Inframe Insertion" = "#380099", "Splice Donor" = "#84714F","SNV Multi-hit" = "#242423")
lgd2_cols <- c("Duplication" = "#CC79A7", "Deletion" = "#D8BFD8", "Inversion" = "#56B4E9", "Translocation" = "#E69F00", "SV Multi-hit" = "#242423")
lgd3_cols <- c( "Duplication" = "#009E73", "Deletion" = "#0072B2" )
lgd4_cols <- c( "Deficient" = "#242423", "Proficient" = "white" )
  
lgd1 <- Legend(at = names(lgd1_cols), legend_gp = gpar(filll = lgd1_cols), title = "SNV")
lgd2 <- Legend(at = names(lgd2_cols), legend_gp = gpar(col = lgd2_cols, lwd = 4, lineend="butt"), type = "lines", title = "SV")
lgd3 <- Legend(at = names(lgd3_cols), legend_gp = gpar(col = lgd3_cols, lwd = 4, lineend="butt"), type = "lines", title = "CNV")

pd <- packLegend(list = list(lgd1, lgd2, lgd3), direction = "vertical")

oncoplot <- oncoPrint(matrix,
                      alter_fun = alter_fun, 
                      col = col, 
                      row_names_side = "left", 
                      show_pct = F,
                      # pct_side = "right", 
                      pct_gp = gpar(fontsize = 16), 
                      row_split = subplots$sub_plot, 
                      gap = unit(5, "mm"),
                      top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(bar_width = 0.8),
                                                         height = unit(2, units = "cm")),
                      right_annotation = rowAnnotation(pathogenic_pc = anno_text(paste0(pathogenic$percentage, "%"), gp = gpar(fontsize = 16)),
                                                       overall_pc = anno_text(paste0("(", overall$percentage, "%)"), gp = gpar(fontsize = 16)),
                                                       rbar = anno_oncoprint_barplot(bar_width = 0.8), 
                                                       width = unit(5.5, units = "cm")))

pdf(file = "panelANEW.pdf", width = 18, height = 6.5)
draw(oncoplot, ht_gap = unit(7, "mm"), row_km = 2, show_heatmap_legend = FALSE, annotation_legend_list = pd)
dev.off()

###########################################
# Driver Gene Mutations Barplot - Panel B #
###########################################
driver_muts_all <- filter(muts, GENE %in% drivers)

barplot <- driver_muts_all %>% 
  mutate(CATEGORY = ifelse(CATEGORY == "Structural Rearrangement", CLUSTER, CATEGORY)) %>%
  dplyr::select(SAMPLE, CATEGORY) %>% distinct()

table <- barplot %>%
  group_by(CATEGORY) %>% tally() %>%
  add_row(CATEGORY = "All", n = length(unique(barplot$SAMPLE))) %>%
  mutate(percentage = (n / 324) * 100,
         CATEGORY = factor(CATEGORY, levels = rev(c("Coding SNV/INDELs", "Simple SV", "Clustered SV", "CNV_DUP", "CNV_DEL", "All"))))

cols <- c("CNV_DEL" = "#0072B2",
          "CNV_DUP" = "#009E73",
          "Coding SNV/INDELs" = "#CC79A7",
          "Clustered SV" = "#E69F00",
          "Simple SV" = "#56B4E9",
          "All" = "#242423")

p1 <- ggplot(table, aes(y = CATEGORY, x = percentage, fill = CATEGORY)) +
  geom_bar(stat = "identity", show.legend = F) +
  xlim(0,100) +
  xlab("Patients with candidte drivers (%)") +
  geom_text(aes(label = round(percentage)),col="white", hjust = 2) +
  scale_y_discrete(labels = c("CNV_DEL" = "CNV deletion", "CNV_DUP" = "CNV duplication", "Coding SNV/INDELs" = "Coding SNV/INDEL",
                            "Clustered SV" = "Clustered SV", "Simple SV" = "Simple SV", "All" = "Any"))+
  scale_fill_manual(values = cols)+
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 13),
        axis.text.y = element_text(size = 13))

category_average <- driver_muts_all %>%
  mutate(CATEGORY = ifelse(CATEGORY == "Structural Rearrangement", CLUSTER, CATEGORY)) %>%
  group_by(SAMPLE, CATEGORY) %>% tally() %>%
  group_by(CATEGORY) %>%
  dplyr::summarise(CATEGORY = unique(CATEGORY),
                   MEAN = mean(n), 
                   SD = sd(n))

all_average <- driver_muts_all %>%
  dplyr::select(SAMPLE, GENE) %>%
  distinct() %>%
  group_by(SAMPLE) %>%
  tally() %>%
  dplyr::summarise(CATEGORY = "All",
                   MEAN = mean(n),
                   SD = sd(n)) 

average <- rbind(category_average, all_average)
average$CATEGORY<- factor(average$CATEGORY, levels = rev(c("Coding SNV/INDELs", "Simple SV", "Clustered SV", "CNV_DUP", "CNV_DEL", "All")))

p2 <- ggplot(average, aes(x = MEAN, y = CATEGORY, color = CATEGORY)) + 
  geom_pointrange(aes(xmin = MEAN-SD, xmax = MEAN+SD), show.legend = F) +
  xlab("Number of candidate drivers") +
  theme_classic() +
  geom_text(aes(label = round(MEAN, 1)), vjust = -1, show.legend = F) +
  scale_color_manual(values = cols) +
  scale_y_discrete(labels=c("CNV_DEL"="CNV deletion","CNV_DUP" = "CNV duplication","Coding SNV/INDELs"="Coding SNV/INDEL",
                            "Clustered SV"="Clustered SV","Simple SV"="Simple SV","All"="Any"))+
  theme(axis.text.y = element_blank(),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_blank())

## Combine plots
pdf(file = "panelBNEW.pdf", width = 9, height = 4)
plot_grid(p1, p2, nrow = 1, align = "hv")
dev.off()

#############################################
# Driver Gene Mutation Categories - Panel C #
#############################################
aggregated <- driver_muts_all %>%
  mutate(CATEGORY = ifelse(CATEGORY == "Structural Rearrangement", CLUSTER, CATEGORY)) %>%
  dplyr::select(SAMPLE, GENE, CATEGORY) %>%
  distinct() %>%
  aggregate(data = ., CATEGORY~., FUN = paste, collapse = ";") %>%
  mutate(CATEGORY = ifelse(grepl(";", CATEGORY), "Multi-Hit", CATEGORY))

sample_tally <- aggregated %>% 
  dplyr::group_by(GENE, CATEGORY) %>% 
  tally()

totals <- aggregated %>%
  dplyr::group_by(GENE) %>%
  dplyr::summarize(total = n(),
                   percentage = paste0(round((total / 324) * 100, 0), "%"))

cols <- c("CNV_DEL" = "#0072B2",
          "CNV_DUP" = "#009E73",
          "Coding SNV/INDELs" = "#CC79A7",
          "Clustered SV" = "#E69F00",
          "Simple SV" = "#56B4E9",
          "Multi-Hit" = "black")

pdf(file = "Fig3C.pdf", width = 10, height = 4)
ggplot(sample_tally, aes(x = reorder(GENE, -n, sum), y = n, fill = CATEGORY)) +
  geom_bar(stat="identity") +
  geom_text(data = totals, aes(x = GENE, total + 30, label = percentage, fill = NULL), angle = 45, size = 5,fontface=2) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size=14, colour="black"),
        axis.text.y = element_text(size=14, colour="black"),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14,face=2)) +
  scale_fill_manual(values = cols,labels=c("CNV_DEL"="CNA deletion","CNV_DUP" = "CNA duplication","Coding SNV/INDELs"="Coding SNV/INDEL",
                            "Clustered SV"="Clustered SV","Simple SV"="Simple SV","Multi-Hit"="Multi-hit"),name="Category")+
  xlab("Candidate driver gene") +
  ylab("Number of patients") +
  scale_x_discrete(expand = expansion(mult = 0.05))
dev.off()
