##########################################
# Script to plot SELECT analysis heatmap #
##########################################

library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)

## Load in alpi file output from select
alpi <- fread("data/selectOutput.txt.gz")

to_plot <- c("CEP89_SV", "CCNE1_SV", "PTEN_SV",
             "TP53_SNV", "NF1_SNV", "RB1_SNV", "SLC35G5_SNV", "BRCA1_SNV", "CDK12_SNV", "TAS2R43_SNV", "BRCA2_SNV",
             "AXIN1_CNV", "PCM1_CNV", "WRN_CNV", "CEP43_CNV", "LEPROTL1_CNV", "DNMT3A_CNV", "ARID1B_CNV", "TCEA1_CNV",
             "HRD_COMPLEX", "ecDNA_COMPLEX", "Chromothripsis_COMPLEX", "WholeGenomeDuplication_COMPLEX",
             "complexI_SNV", "complexIII_SNV", "complexIV_SNV", "complexV_SNV")

## Filter alpi file to only the genes you are interested in 
alpi <- filter(alpi, SFE_1 %in% to_plot & SFE_2 %in% to_plot)

## Add in co-occurrence / mutual-exclusivity results
alpi$result <- ifelse(alpi$APC >= 0.003 & alpi$direction == "ME", "ME",
                           ifelse(alpi$APC >= 0.003 & alpi$direction == "CO", "CO",
                                  ifelse(alpi$APC < 0.003 & alpi$direction == "ME", "W_ME",
                                         ifelse(alpi$APC < 0.003 & alpi$direction == "CO", "W_CO", NA))))

## Format results so we have data for all possible combinations
matrix <- alpi %>%
  dplyr::select(SFE_1, SFE_2, result) %>%
  add_row(SFE_1 = setdiff(to_plot, alpi$SFE_1)) %>% 
  add_row(SFE_2 = setdiff(to_plot, alpi$SFE_2)) %>%
  complete(SFE_1, SFE_2, explicit = F) %>% 
  filter(!is.na(SFE_1) & !is.na(SFE_2))

matrix$result <- apply(matrix, 1, function(x) { ifelse(is.na(x[3]), 
                                                       purrr::discard(matrix$result[(matrix$SFE_1 == x[1] & matrix$SFE_2 == x[2]) | (matrix$SFE_1 == x[2] & matrix$SFE_2 == x[1])], is.na),
                                                       x[3]) })

matrix$result <- ifelse(matrix$SFE_1 == matrix$SFE_2, "SELF", matrix$result)
matrix$result[is.na(matrix$result)] <- "NONE"

matrix$int_type <- apply(matrix, 1, function(x) { paste0(strsplit(x[1], "_")[[1]][2], "-", strsplit(x[2], "_")[[1]][2]) })

## If the comparison is to itself, find the number of events for this gene
events <- unique(rbind(dplyr::select(alpi, SFE_1, support_1), dplyr::select(alpi, SFE_2, support_2), use.names = F))
colnames(events) <- c("gene", "count")
matrix$cluster <- apply(matrix, 1, function(x) { ifelse(x[3] == "SELF", events$count[events$gene == x[1]], "") })

## If the comparison has a direction, find the overlap and predicted overlap
alpi$combined_overlap <- paste0(alpi$overlap, "\n", round(alpi$r_overlap, 1))
matrix$overlap <- apply(matrix, 1, function(x) { ifelse(x[3] != "NONE" & x[3] != "SELF", unique(alpi$combined_overlap[(alpi$SFE_1 == x[1] & alpi$SFE_2 == x[2]) | (alpi$SFE_1 == x[2] & alpi$SFE_2 == x[1])]), "") })

## Put genes / events into groups
set1 <- sort(c("complexI_SNV", "complexIII_SNV", "complexIV_SNV", "complexV_SNV"))
set2 <- sort(c("TP53_SNV", "NF1_SNV", "RB1_SNV", "SLC35G5_SNV", "BRCA1_SNV", "CDK12_SNV", "TAS2R43_SNV", "BRCA2_SNV"))
set3 <- sort(c("CEP89_SV", "CCNE1_SV", "PTEN_SV"))
set4 <- sort(c("AXIN1_CNV", "PCM1_CNV", "WRN_CNV", "CEP43_CNV", "LEPROTL1_CNV", "DNMT3A_CNV", "ARID1B_CNV", "TCEA1_CNV"))
set5 <- sort(c("HRD_COMPLEX", "ecDNA_COMPLEX", "Chromothripsis_COMPLEX", "WholeGenomeDuplication_COMPLEX"))

levels <- c(set1, set2, set3, set4, set5)

ordering <- levels %>% as.data.frame() %>% dplyr::mutate(order = 1:n())
matrix <- left_join(matrix, ordering, by = c("SFE_1" = "."))
matrix <- left_join(matrix, ordering, by = c("SFE_2" = "."))
matrix <- matrix %>% arrange(order.x, order.y)

matrix <- matrix %>% mutate(SFE_1 = factor(SFE_1, levels = levels),
                            SFE_2 = factor(SFE_2, levels = rev(levels)))

num_blank <- 0
labels <- c()
for(i in 1:length(unique(matrix$SFE_1))){
  add <- c(rep("KEEP", num_blank), rep("REMOVE", length(unique(matrix$SFE_1))-num_blank))
  labels <- c(labels, add)
  num_blank <- num_blank + 1
}

matrix$labels <- labels
matrix$overlap <- ifelse(matrix$labels == "REMOVE", "", matrix$overlap)

matrix$label_colour <- ifelse(matrix$result == "CO" | matrix$result == "ME", "light", "dark")

## Plot the heatmap
pdf("SELECT.pdf", width = 10, height = 10)
ggplot(matrix, aes(SFE_1, SFE_2)) +
  geom_tile(aes(fill = result), colour = "#F2F2F2", size = 0.5) +
  theme(panel.background = element_blank(),
        axis.text.x = element_text(angle = -45, hjust=1, size = 7),
        axis.text.y = element_text(size = 7),
        plot.margin = unit(c(1, 1, 2, 1), "cm"),
        legend.title = element_blank()) +
  scale_x_discrete(position = "top") +
  coord_equal() + 
  xlab("") +
  ylab("") +
  geom_text(aes(label = cluster), size = 2.5, color = "#666666", fontface = "bold") +
  geom_text(aes(label = overlap, color = label_colour), size = 1.5, show.legend = F, fontface = "bold") +
  scale_color_manual(values = c("dark" = "#201C1C", "light" = "white")) +
  scale_fill_manual(labels = c("Mutually exclusive", "Co-occurring", "Weakly mutually-exclusive", "Weakly co-occurring","No pattern", "Not tested"),
                    values = c("ME" = "#892BE1", "CO" = "#008002", 
                               "W_ME" = "#DCBFF7", "W_CO" = "#BFF3BF", 
                               "NONE" = "white", "SELF" = "grey90"))
dev.off()
