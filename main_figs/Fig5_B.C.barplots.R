###################################
# Mitochondria Figure - Panel B/C #
###################################

library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)


###############################################################
# Mitochondrial Respiratory Complexes Mutational Consequences #
###############################################################
complexes <- fread("data/mtComplexes.csv")

consequence_dict <- c("frameshift_variant" = "Frameshift Variant",
                      "stop_gained" = "Stop Gained",
                      "missense_variant" = "Missense Variant",
                      "stop_retained" = "Stop Retained",
                      "inframe_deletion" = "Inframe Deletion", 
                      "inframe_insertion" = "Inframe Insertion",
                      "synonymous_variant" = "Synonymous Variant")

complex_dict <- c("complexI" = "Complex I",
                  "complexIII" = "Complex III",
                  "complexIV" = "Complex IV",
                  "complexV" = "Complex V")

muts <- fread("data/allCallsFINALNEW.txt.gz") %>%
  filter(somatic_status == "Somatic") %>%
  separate_rows(all_of(names(.)[which(grepl(";", .))]), sep = ";", convert = T) %>%
  merge(., complexes, by.x = "SYMBOL", by.y = "gene")  %>%
  mutate(Consequence = str_replace_all(Consequence, pattern = consequence_dict),
         complex = str_replace_all(complex, pattern = complex_dict))

cols <- c("Frameshift Variant" = "#0B6FE3",
          "Stop Gained" = "#FBAC11",
          "Missense Variant" = "#137547",
          "Stop Retained" = "#C84630",
          "Inframe Deletion" = "#FFC6AC", 
          "Inframe Insertion" = "#A1CCA5",
          "Synonymous Variant" = "#7C77B9")

pdf("complex_mutational_consequences.pdf", width = 8.5, height = 4)
ggplot(muts, aes(x = complex, fill = Consequence)) +
  geom_bar() +
  xlab("Mitochondrial Respiratory Chain Complex") +
  ylab("Number of Mutations") +
  coord_cartesian(ylim = c(0,450)) +
  theme(#legend.position = "none",
    panel.border = element_rect(colour = "#242423", fill = NA, size = rel(0.7)),
    axis.ticks = element_line(colour = "#242423",  size = rel(0.7), lineend = "round"),
    axis.text = element_text(colour = "#242423", size = rel(1.1)),
    axis.title = element_text(colour = "#242423", size = rel(1)),
    title = element_text(colour = "#242423", size = rel(1.05)),
    strip.background = element_blank(),
    axis.line = element_blank(),
    panel.grid.major = element_line(colour = "gray95"),
    panel.grid.minor = element_line(colour = "gray95"), 
    panel.background = element_blank(), 
    legend.title = element_text(face = "bold")) +
  scale_fill_manual(values = cols)
dev.off()

###############################################
# Mitochondrial Genes Mutational Consequences #
###############################################
consequence_dict <- c("frameshift_variant" = "Frameshift Variant",
                      "stop_gained" = "Stop Gained",
                      "missense_variant" = "Missense Variant",
                      "stop_retained" = "Stop Retained",
                      "inframe_deletion" = "Inframe Deletion", 
                      "inframe_insertion" = "Inframe Insertion",
                      "synonymous_variant" = "Synonymous Variant",
                      "intergenic_variant" = "Intergenic Variant",
                      "non_coding_transcript_exon_variant" = "Non-Coding Transcript Exon Variant")

muts <- fread("data/allCallsFINALNEW.txt.gz") %>%
  filter(somatic_status == "Somatic") %>%
  separate_rows(all_of(names(.)[which(grepl(";", .))]), sep = ";", convert = T) %>%
  mutate(Consequence = str_replace_all(Consequence, pattern = consequence_dict),
         SYMBOL = str_replace_all(SYMBOL, pattern = c("^-$" = "NCR"))) %>%
  mutate(SYMBOL = factor(SYMBOL, levels = c("NCR", sort(unique(SYMBOL[SYMBOL != "NCR"])))))

cols <- c("Frameshift Variant" = "#0B6FE3",
          "Stop Gained" = "#FBAC11",
          "Missense Variant" = "#137547",
          "Stop Retained" = "#C84630",
          "Inframe Deletion" = "#FFC6AC", 
          "Inframe Insertion" = "#A1CCA5",
          "Synonymous Variant" = "#7C77B9",
          "Intergenic Variant" = "#646E78",
          "Non-Coding Transcript Exon Variant" = "#2297E6")

pdf("genes_mutational_consequences.pdf", width = 11, height = 4)
ggplot(muts, aes(x = SYMBOL, fill = Consequence)) +
  geom_bar() +
  xlab("Mitochondrial Feature") +
  ylab("Number of Mutations") +
  coord_cartesian(ylim = c(0,155)) +
  theme(#legend.position = "none",
    panel.border = element_rect(colour = "#242423", fill = NA, size = rel(0.7)),
    axis.ticks = element_line(colour = "#242423",  size = rel(0.7), lineend = "round"),
    axis.text = element_text(colour = "#242423", size = rel(1.1)),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(colour = "#242423", size = rel(1)),
    title = element_text(colour = "#242423", size = rel(1.05)),
    strip.background = element_blank(),
    axis.line = element_blank(),
    panel.grid.major = element_line(colour = "gray95"),
    panel.grid.minor = element_line(colour = "gray95"), 
    panel.background = element_blank(),
    legend.title = element_text(face = "bold")) +
  scale_fill_manual(values = cols) +
  scale_x_discrete(expand = expansion(mult = 0.03))
dev.off()
