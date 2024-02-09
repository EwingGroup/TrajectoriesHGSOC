##################################################################################
# Figure showing the large number of SNVs/SVs/CNVs hitting genes in each patient #
##################################################################################

library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)

long <- fread("data/longFormatPathogenicOnly.txt.gz")

genes <- fread("data/ensemblGenes107.bed")
genes <- filter(genes, V5 == "protein_coding")
genes_number <- nrow(genes)

##########
# Plot 1 #
##########
percentage_of_protein_coding <- long %>% 
  mutate(CATEGORY = ifelse(CATEGORY == "CNV", MUT_TYPE, CATEGORY)) %>% 
  dplyr::select(GENE, CATEGORY) %>% 
  distinct() %>% 
  group_by(CATEGORY) %>% 
  tally() %>% 
  mutate(percentage = (n / genes_number) * 100)

percentage_of_protein_coding_overall <- long %>% 
  dplyr::select(GENE) %>% 
  distinct() %>% 
  tally() %>% 
  mutate(CATEGORY = "All", percentage = (n / genes_number) * 100) %>% 
  dplyr::select(CATEGORY, n, percentage)

all_categories <- rbind(percentage_of_protein_coding, percentage_of_protein_coding_overall)

average_per_sample <- long %>% 
  mutate(CATEGORY = ifelse(CATEGORY == "CNV", MUT_TYPE, CATEGORY)) %>%
  dplyr::select(SAMPLE, GENE, CATEGORY) %>% distinct() %>%
  group_by(SAMPLE, CATEGORY) %>%
  tally() %>% 
  group_by(CATEGORY) %>% 
  dplyr::summarise(CATEGORY = unique(CATEGORY), MEAN = mean(n), SD = sd(n))

average_per_sample_overall <- long %>% 
  dplyr::select(SAMPLE, GENE) %>% distinct() %>%
  group_by(SAMPLE) %>% 
  tally() %>%
  dplyr::summarise(CATEGORY = "All", MEAN = mean(n), SD = sd(n))

all_averages <- rbind(average_per_sample, average_per_sample_overall)

all_categories$CATEGORY <- factor(all_categories$CATEGORY, levels = all_averages[order(-all_averages$MEAN),]$CATEGORY)
all_averages$CATEGORY <- factor(all_averages$CATEGORY, levels = all_averages[order(-all_averages$MEAN),]$CATEGORY)

cols <- c("CNV_DEL" = "#0072B2",
          "CNV_DUP" = "#009E73",
          "SNV" = "#CC79A7",
          "SV" = "#E69F00",
          "All" = "#242423")

p1 <- ggplot(all_categories, aes(x = percentage, y = CATEGORY, fill = CATEGORY)) +
  geom_bar(stat = "identity", show.legend = F) +
  ylab("") +
  xlab("Protein coding genes (%)") +
  xlim(0, 100) +
  scale_fill_manual(values = cols) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 18,colour="black"),
        axis.text.x = element_text(size = 14,colour="black"),
        axis.title.x = element_text(size = 18,colour="black"))+
  scale_y_discrete(labels=c("CNV_DEL" = "CNA deletion", "CNV_DUP" = "CNA duplication","SNV" = "Coding SNV/INDEL",
                    "SV" = "Structural Variant", "All" = "Any")) +
  geom_text(aes(label = round(percentage)),col="white", hjust = 2, size=5, fontface=2)

p2 <- ggplot(all_averages, aes(x = MEAN, y = CATEGORY, color = CATEGORY)) + 
  geom_pointrange(aes(xmin = MEAN-SD, xmax = MEAN+SD), show.legend = F) +
  xlab("Disrupted genes per sample") +
  theme_classic() +
  geom_text(aes(label = round(MEAN, 1)), vjust = -1, show.legend = F,size=5,fontface=2) +
  scale_color_manual(values = cols) +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size=14,colour="black"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 18))

pdf("individual_panels/main_figs/Fig1/Figure1_B.pdf", width = 12, height = 4)
plot_grid(p1, p2, nrow = 1, align = "hv", axis = "lrtb")
dev.off()

##########
# Plot 2 #
##########
top_mutated_categories <- long %>% 
  mutate(CATEGORY = ifelse(CATEGORY == "CNV", MUT_TYPE, CATEGORY)) %>% 
  dplyr::select(SAMPLE, GENE, CATEGORY) %>% distinct() %>% 
  group_by(SAMPLE, GENE) %>%
  dplyr::mutate(CATEGORY = ifelse(length(SAMPLE) > 1, "multi_hit", CATEGORY)) %>%
  distinct() %>% 
  group_by(GENE, CATEGORY) %>% 
  tally() %>% 
  group_by(GENE) %>% 
  dplyr::mutate(total = sum(n))

top_mutated <- top_mutated_categories %>% 
  dplyr::select(GENE, total) %>% 
  distinct() %>% 
  arrange(-total) %>% 
  filter(GENE != "NULL") %>% 
  head(n = 20)

to_plot <- filter(top_mutated_categories, GENE %in% top_mutated$GENE)

cols <- c("CNV_DEL" = "#0072B2",
          "CNV_DUP" = "#009E73",
          "SNV" = "#CC79A7",
          "SV" = "#E69F00",
          "multi_hit" = "#242423")

to_plot$percentage <- paste0(round((to_plot$total / 324) * 100), "%")

pdf("individual_panels/main_figs/Fig1/Figure1_C.pdf", width = 10, height = 4)
ggplot(to_plot, aes(x = reorder(GENE, -total), y = n, fill = CATEGORY)) +
  geom_bar(stat = "identity") +
  xlab("Most frequently disrupted genes") +
  ylab("Number of patients") +
  scale_x_discrete(expand = expansion(mult = 0.05)) +
  theme_classic() +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size=14,colour="black"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size=14, colour="black"),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18)) +
  scale_fill_manual(values = cols, labels=c("CNV_DEL" = "CNA deletion", "CNV_DUP" = "CNA duplication","SNV" = "Coding SNV/INDEL",
                                            "SV" = "Structural Variant", "multi_hit" = "Multi-Hit")) +
  geom_text(data = to_plot, aes(x = GENE, total + 30, label = percentage, fill = NULL), angle = 45, size = 5)
dev.off()
