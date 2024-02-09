#######################################
# Driver Genes Split By Genomic State #
#######################################

library(cowplot)
library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

## percentages should be as a percentage of the genomic state not of 324.



complex <- fread("data/SummaryTableSampleCallOfcSV.tsv", header = T) %>%
  mutate(Sample = gsub("-P$", "", Sample),
         HRD_WGD = paste0(HRDetectStat, " - ", WholeGenomeDuplication))

complex$HRD_WGD[complex$HRD_WGD == "Absent - Absent"] <- "Neither"
complex$HRD_WGD[complex$HRD_WGD == "Absent - Present"] <- "WGD Only"
complex$HRD_WGD[complex$HRD_WGD == "Present - Absent"] <- "HRD Only"
complex$HRD_WGD[complex$HRD_WGD == "Present - Present"] <- "HRD+WGD"

counts <- complex %>%
  group_by(HRD_WGD) %>%
  tally()
counts

sv_drivers <- c("CEP89", "CCNE1", "PTEN")
snv_drivers <- c("TP53", "NF1", "RB1", "SLC35G5", "BRCA1", "CDK12", "TAS2R43", "BRCA2") 
cnv_drivers <- c("AXIN1", "PCM1", "WRN", "CEP43", "LEPROTL1", "DNMT3A", "ARID1B", "TCEA1")
drivers <- unique(c(sv_drivers, snv_drivers, cnv_drivers))

driver_muts_all <- fread("data/oncoplotCombinedExonsHit.txt.gz") %>%
  filter(GENE %in% drivers) %>%
  left_join(complex, by = c("SAMPLE" = "Sample"))

barplot <- driver_muts_all %>% 
  mutate(CATEGORY = ifelse(CATEGORY == "Structural Rearrangement", CLUSTER, CATEGORY)) %>%
  dplyr::select(SAMPLE, HRD_WGD, CATEGORY) %>% distinct()

counts2 <- counts %>% rename("sample_n" = "n")

table <- barplot %>%
  group_by(CATEGORY, HRD_WGD) %>% tally() %>% ungroup() %>%
  left_join(., counts2, by = "HRD_WGD") %>%
  mutate(percentage = (n / sample_n) * 100,
         CATEGORY = factor(CATEGORY, levels = rev(c("Coding SNV/INDELs", "Simple SV", "Clustered SV", "CNV_DUP", "CNV_DEL"))),
         HRD_WGD = factor(HRD_WGD, levels = rev(c("Neither", "WGD Only", "HRD Only", "HRD+WGD"))))

cols <- c("CNV_DEL" = "#0072B2",
          "CNV_DUP" = "#009E73",
          "Coding SNV/INDELs" = "#CC79A7",
          "Clustered SV" = "#E69F00",
          "Simple SV" = "#56B4E9",
          "All" = "#242423")

p1 <- ggplot(table, aes(y = CATEGORY, x = percentage, fill = CATEGORY)) +
  geom_bar(stat = "identity", position = "dodge", show.legend = F) +
  xlim(0,100) +
  xlab("Patients with\ncandidate drivers (%)") +
  facet_wrap(facets = vars(HRD_WGD), ncol = 1, strip.position = "left") +
  geom_text(aes(label = round(percentage)), col="white", hjust = 1.5, position = position_dodge(width = 0.9)) +
  scale_y_discrete(labels = c("CNV_DEL" = "CNV deletion", "CNV_DUP" = "CNV duplication", "Coding SNV/INDELs" = "Coding SNV/INDEL",
                              "Clustered SV" = "Clustered SV", "Simple SV" = "Simple SV", "All" = "Any"), 
                   expand = expansion(add = 0.8))+
  scale_fill_manual(values = cols)+
  theme_classic() +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "#242423", fill = NA, size = rel(0.7)),
        axis.ticks = element_line(colour = "#242423",  size = rel(0.7), lineend = "round"),
        axis.text = element_text(colour = "#242423", size = rel(1.1)),
        axis.title = element_text(colour = "#242423", size = rel(1)),
        title = element_text(colour = "#242423", size = rel(1.05)), 
        strip.text = element_text(colour = "#242423", size = rel(1.05), face = "bold"),
        strip.background = element_rect(color = "black", fill = "gray95"),
        axis.line = element_blank(),
        panel.grid.major = element_line(colour = "gray95"),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.title.y = element_blank())

average <- driver_muts_all %>%
  mutate(CATEGORY = ifelse(CATEGORY == "Structural Rearrangement", CLUSTER, CATEGORY)) %>%
  group_by(SAMPLE, HRD_WGD, CATEGORY) %>% tally() %>%
  group_by(CATEGORY, HRD_WGD) %>%
  dplyr::summarise(CATEGORY = unique(CATEGORY),
                   HRD_WGD = unique(HRD_WGD),
                   MEAN = mean(n), 
                   SD = sd(n), .groups = "drop") %>%
  mutate(CATEGORY = factor(CATEGORY, levels = rev(c("Coding SNV/INDELs", "Simple SV", "Clustered SV", "CNV_DUP", "CNV_DEL"))),
         HRD_WGD = factor(HRD_WGD, levels = rev(c("Neither", "WGD Only", "HRD Only", "HRD+WGD"))))

p2 <- ggplot(average, aes(x = MEAN, y = CATEGORY, color = CATEGORY)) + 
  geom_pointrange(aes(xmin = MEAN-SD, xmax = MEAN+SD), show.legend = F, position = position_dodge(width = 0.9)) +
  xlab("Number of\ncandidate drivers") +
  facet_wrap(facets = vars(HRD_WGD), ncol = 1, strip.position = "right") +
  theme_classic() +
  geom_text(aes(x = MEAN+SD, label = round(MEAN, 1)), show.legend = F, position = position_dodge(width = 0.9), hjust = -0.5) +
  scale_color_manual(values = cols) +
  scale_x_continuous(limits=c(-10,40)) +
  scale_y_discrete(labels=c("CNV_DEL"="CNV deletion","CNV_DUP" = "CNV duplication","Coding SNV/INDELs"="Coding SNV/INDEL",
                            "Clustered SV"="Clustered SV","Simple SV"="Simple SV","All"="Any"), 
                   expand = expansion(add = 0.8))+
  theme(legend.position = "none", 
        panel.border = element_rect(colour = "#242423", fill = NA, size = rel(0.7)),
        axis.ticks = element_line(colour = "#242423",  size = rel(0.7), lineend = "round"),
        axis.text = element_text(colour = "#242423", size = rel(1.1)),
        axis.title = element_text(colour = "#242423", size = rel(1)),
        title = element_text(colour = "#242423", size = rel(1.05)),
        strip.background.y = element_blank(),
        strip.text = element_blank(),
        axis.line = element_blank(),
        panel.grid.major = element_line(colour = "gray95"),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())

# ## Combine plots
# pdf(file = "~/Desktop/panelBNEW.pdf", width = 6, height = 6)
# plot_grid(p1, p2, nrow = 1, align = "h", rel_widths = c(0.62,0.38))
# dev.off()

#############################################
# Driver Gene Mutation Categories - Panel C #
#############################################
aggregated <- driver_muts_all %>%
  mutate(CATEGORY = ifelse(CATEGORY == "Structural Rearrangement", CLUSTER, CATEGORY)) %>%
  dplyr::select(SAMPLE, HRD_WGD, GENE, CATEGORY) %>%
  distinct() %>%
  aggregate(data = ., CATEGORY~., FUN = paste, collapse = ";") %>%
  mutate(CATEGORY = ifelse(grepl(";", CATEGORY), "Multi-Hit", CATEGORY))

sample_tally <- aggregated %>% 
  dplyr::group_by(GENE, HRD_WGD, CATEGORY) %>% 
  tally()

totals <- aggregated %>%
  dplyr::group_by(GENE, HRD_WGD) %>%
  dplyr::summarize(total = n(), .groups = "drop") %>%
  left_join(., counts, by = "HRD_WGD") %>%
  mutate(percentage = paste0(round((total / n) * 100, 0), "%"))

p3 <- ggplot(sample_tally, aes(x = GENE, y = n, fill = CATEGORY)) +
  geom_bar(stat="identity") +
  facet_wrap(facets = vars(HRD_WGD), ncol = 1, strip.position = "right") +
  scale_x_discrete(expand = expansion(add = 0.8)) +
  geom_text(data = totals, aes(x = GENE, total + 30, label = percentage, fill = NULL), angle = 45, size = 3) +
  labs(x = "Candidate driver gene", y = "Number of Patients") +
  coord_cartesian(ylim = c(0,175)) +
  scale_fill_manual(values = cols) +
  theme_classic() +
  theme(legend.position = "none", 
        panel.border = element_rect(colour = "#242423", fill = NA, size = rel(0.7)),
        axis.ticks = element_line(colour = "#242423",  size = rel(0.7), lineend = "round"),
        axis.text = element_text(colour = "#242423", size = rel(1.1)),
        axis.title = element_text(colour = "#242423", size = rel(1)),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        title = element_text(colour = "#242423", size = rel(1.05)),
        strip.background.y = element_blank(),
        strip.text = element_blank(),
        axis.line = element_blank(),
        panel.grid.major = element_line(colour = "gray95"),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank())

pdf(file = "~/Desktop/supFig11.pdf", width = 11, height = 6)
plot_grid(p1, p2, p3, nrow = 1, align = "h", rel_widths = c(3,1.5,5))
dev.off()

  
  
  
  
  
  

  
