##################################################
# Boxplots Comparing Elements of SELECT Analysis #
##################################################

library(data.table)
library(dplyr)
library(ggplot2)
library(rstatix)
library(ggpubr)



muts <- fread("data/allCallsFINALNEW.txt.gz") %>%
  filter(somatic_status == "Somatic") %>%
  mutate(sample = gsub("_2$|-P$", "", sample)) %>%
  group_by(sample) %>% tally()

complex <- fread("data/SummaryTableSampleCallOfcSV.tsv", header = T) %>%
  dplyr::select(-V1, -Cohort, -dplyr::contains("Starfish")) %>%
  mutate(Sample = gsub("-P$", "", Sample), SS_Chromothripsis = ifelse(SS_Chromothripsis != "Absent", "Present", SS_Chromothripsis)) %>%
  dplyr::rename(sample = "Sample") 

combined <- left_join(complex, muts, by = "sample")
combined$n[is.na(combined$n)] <- 0 
combined$HRD_WGD <- paste0(combined$HRDetectStat, " - ", combined$WholeGenomeDuplication)

combined$HRD_WGD[combined$HRD_WGD == "Absent - Absent"] <- "Neither"
combined$HRD_WGD[combined$HRD_WGD == "Absent - Present"] <- "WGD Only"
combined$HRD_WGD[combined$HRD_WGD == "Present - Absent"] <- "HRD Only"
combined$HRD_WGD[combined$HRD_WGD == "Present - Present"] <- "HRD+WGD"

pdf("WGD_mtMuts.pdf", width = 6, height = 4)
ggplot(combined, aes(x = WholeGenomeDuplication, y = n)) +
  geom_jitter(width = 0.2) +
  geom_boxplot(outlier.shape = NA, fill = NA) +
  geom_signif(comparisons = list(c(1, 2)), y_position = 31, textsize = 3, size = 0.25, tip_length = NA, map_signif_level = F, test = "wilcox.test") +
  coord_cartesian(ylim = c(0,35)) +
  labs(x = "Whole Genome Duplication", y = "Number of Mitochondrial Somatic Mutations") +
  theme_classic()
dev.off()

pdf("HRD_mtMuts.pdf", width = 6, height = 4)
ggplot(combined, aes(x = HRDetectStat, y = n)) +
  geom_jitter(width = 0.2) +
  geom_boxplot(outlier.shape = NA, fill = NA) +
  geom_signif(comparisons = list(c(1, 2)), y_position = 31, textsize = 3, size = 0.25, tip_length = NA, map_signif_level = F, test = "wilcox.test") +
  coord_cartesian(ylim = c(0,35)) +
  labs(x = "Homologous Recombination Deficiency", y = "Number of Mitochondrial Somatic Mutations") +
  theme_classic()
dev.off()

stat.test <- combined %>% 
  wilcox_test(n ~ HRD_WGD) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "HRD_WGD")

pdf("HRD_WGD_mtMuts.pdf", width = 6, height = 4)

combined$HRD_WGD<- factor(combined$HRD_WGD,levels=c("HRD+WGD","HRD Only","WGD Only","Neither"))
stat.test <- combined %>% 
  wilcox_test(n ~ HRD_WGD) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "HRD_WGD")
ggplot(combined, aes(x = HRD_WGD, y = n)) +
  geom_jitter(width = 0.2) +
  geom_boxplot(outlier.shape = NA, fill = NA) +
 # stat_pvalue_manual(stat.test, label = "p = {p}", size = 3, label.size = 3, tip.length = 0) +
 # coord_cartesian(ylim = c(0,40)) +
  labs(x = "Genomic state", y = "Number of Mitochondrial Somatic Mutations") +
  theme_classic()
dev.off()









