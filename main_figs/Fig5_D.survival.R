##############################################################################
# Basic HGSOC survival plot splitting samples into deleterious and wild-type #
##############################################################################

library(ggplot2)
library(ggfortify)
library(survival)
library(dplyr)
library(data.table)
library(tidyr)


muts <- fread("data/allCallsFINALNEW.txt.gz") %>%
  filter(somatic_status == "Somatic") %>%
  separate_rows(all_of(names(.)[which(grepl(";", .))]), sep = ";", convert = T) %>%
  mutate(SIFT = gsub("\\(.*\\)", "", SIFT), PolyPhen = gsub("\\(.*\\)", "", PolyPhen),
         sample = gsub("_2$|-P$", "", sample)) %>%
  filter(IMPACT == "HIGH" | (SIFT == "deleterious" & PolyPhen == "probably_damaging")) %>%
  filter(Consequence != "inframe_deletion" | Consequence != "inframe_insertion")

hrd <- fread("data/HRDetect_results_BRCAstatus.txt") %>%
  mutate(HRDeficiency = ifelse(HRDetect > 0.7, "Deficient", "Proficient"),
         Sample = gsub("_2$|-P$", "", Sample))

clinical <- fread("data/HGSOC_cohorts_minimal_clin.txt") %>%
  mutate(Sample = gsub("_2$|-P$", "", Sample),
         donor_survival_time = ifelse(donor_survival_time > 3650, 3650, donor_survival_time),
         os_event = ifelse(donor_survival_time > 3650, 0, os_event),
         mutation = ifelse(Sample %in% muts$sample, "B_MUTATED", "A_WILDTYPE")) %>%
  left_join(., hrd, by = "Sample")

res.cox <- coxph(Surv(donor_survival_time, os_event) ~ mutation + HRDeficiency + donor_age_at_diagnosis + factor(stage) + strata(cohort), data = clinical)
summary(res.cox) 

km_fit <- survfit(Surv(donor_survival_time, os_event) ~ mutation, data = clinical)
km_fit

pdf("hrd_adjusted.pdf", width = 5, height = 5.25/2)
autoplot(km_fit, surv.size = 1, conf.int.alpha = 0.12) +
  labs(x = "Time (days)", y = "Overall Survival") + theme_bw() +
  coord_cartesian(ylim = c(0,1), xlim = c(0,3650)) +
  scale_color_manual(values = c("A_WILDTYPE" = "#2176FF", "B_MUTATED" = "#FC440F")) +
  scale_fill_manual(values = c("A_WILDTYPE" = "#2176FF", "B_MUTATED" = "#FC440F")) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "#242423", size = rel(0.7)),
        axis.ticks = element_line(colour = "#242423",  size = rel(0.7), lineend = "round"),
        axis.text = element_text(colour = "#242423", size = rel(1.1)),
        axis.title = element_text(colour = "#242423", size = rel(1)),
        title = element_text(colour = "#242423", size = rel(1.05))) 
dev.off()
