library(dplyr)
library(ggplot2)
# library(ggpubr)
library(multcomp)
library(tidyr)
library(yardstick)
library(tidyr)
library(pROC)
library(wrappedtools)
library(stringr)
library(see)
library(ggrepel)
source("myfuns.R")

####### ======  Read and consolidate results data ======== #######

# Thresholds
intthresESM2 <- read.csv("../results/esm2/best_mcc_thresholds.csv", header = TRUE) %>%
  rename(Pathogen  = Dataset, 
         Predictor = Method) %>%
  mutate(Embedder= "esm2")

intthresESM1 <- read.csv("../results/esm1b/best_mcc_thresholds.csv", header = TRUE) %>%
  rename(Pathogen  = Dataset, 
         Predictor = Method) %>%
  mutate(Embedder= "esm1b")

extpreds <- c("bepipred3", "epidope", "epitopevec", "epitope1d")
extthres <- data.frame(Predictor = rep(extpreds, 
                                       each = length(unique(intthresESM2$Pathogen))),
                       Pathogen  = rep(unique(intthresESM2$Pathogen), 
                                       times = length(extpreds)),
                       Threshold = rep(c(0.1512, 0.818, 0.5, 0.5), 
                                       each = length(unique(intthresESM2$Pathogen))),
                       Embedder = "other")

pthres <- bind_rows(intthresESM1, intthresESM2, extthres) %>%
  mutate(Pathogen = gsub("\ ", "\\_", tolower(Pathogen)),
         Predictor = gsub("-", "", tolower(Predictor)),
         Predictor = gsub("esm1b|esm2", "baseline", Predictor))


# Test/Holdout fold information
testFolds <- read.csv("../data/TrainTestFolds.csv", header = TRUE) %>%
  mutate(Pathogen = gsub("\ ", "\\_", tolower(Pathogen)))

# Consolidated results
X <- readRDS("../data/results/consolidated_results_v9.rds") %>%
  filter(!is.na(Class)) %>%
  left_join(pthres, by = c("Predictor", "Pathogen", "Embedder")) %>%
  left_join(testFolds, by = c("Pathogen", "Fold")) %>%
  mutate(Pred = factor(Prob >= Threshold, 
                       levels = c("TRUE", "FALSE"), 
                       labels = c("TRUE", "FALSE"),
                       ordered = TRUE),
         Class = factor(Class == 1, 
                        levels = c("TRUE", "FALSE"), 
                        labels = c("TRUE", "FALSE"),
                        ordered = TRUE),
         Predictor_full = ifelse(Embedder == "other", Predictor, paste0(Predictor, "_", Embedder))) %>%
  dplyr::select(-Threshold) %>%
  as_tibble()




####### ======  Calculate performance metrics ======== #######

Xmetrics.HO <- X %>%
  filter(Holdout == 1) %>%
  dplyr::select(-Holdout, -Fold) %>%
  group_by(Predictor_full, Pathogen) %>%
  arrange(Class) %>%
  summarise(AUC = as.numeric(pROC::roc(response = Class, 
                                       predictor = Prob, 
                                       levels = c("FALSE", "TRUE"))$auc),
            F1  = f_meas_vec(truth = Class, estimate = Pred, event_level = "first"),
            MCC = mcc_vec(truth = Class, estimate = Pred, event_level = "first"), 
            PPV = ppv_vec(truth = Class, estimate = Pred, event_level = "first"),
            NPV = npv_vec(truth = Class, estimate = Pred, event_level = "first"),
            Sensitivity = sens_vec(truth = Class, estimate = Pred, event_level = "first"),
            Specificity = spec_vec(truth = Class, estimate = Pred, event_level = "first"),
            `Balanced accuracy` = (Sensitivity + Specificity) / 2,
            Accuracy = accuracy_vec(truth = Class, estimate = Pred, event_level = "first"),
            .groups = "drop") %>%
  pivot_longer(AUC:Accuracy, names_to = "Metric", values_to = "Value") %>%
  mutate(Set = "Holdout")

Xmetrics.noLeak <- X %>%
  filter(InTraining == 0) %>%
  dplyr::select(-Holdout, -Fold) %>%
  group_by(Predictor_full, Pathogen) %>%
  arrange(Class) %>%
  summarise(AUC = as.numeric(pROC::roc(response = Class, 
                                       predictor = Prob, 
                                       levels = c("FALSE", "TRUE"))$auc),
            F1  = f_meas_vec(truth = Class, estimate = Pred, event_level = "first"),
            MCC = mcc_vec(truth = Class, estimate = Pred, event_level = "first"), 
            PPV = ppv_vec(truth = Class, estimate = Pred, event_level = "first"),
            NPV = npv_vec(truth = Class, estimate = Pred, event_level = "first"),
            Sensitivity = sens_vec(truth = Class, estimate = Pred, event_level = "first"),
            Specificity = spec_vec(truth = Class, estimate = Pred, event_level = "first"),
            `Balanced accuracy` = (Sensitivity + Specificity) / 2,
            Accuracy = accuracy_vec(truth = Class, estimate = Pred, event_level = "first"),
            .groups = "drop") %>%
  pivot_longer(AUC:Accuracy, names_to = "Metric", values_to = "Value") %>%
  mutate(Set = "NoLeakage")

Xmetrics <- bind_rows(Xmetrics.HO, Xmetrics.noLeak) %>%
  rename(Dataset = Pathogen, Method = Predictor_full) %>%
  mutate(Dataset = fix_dataset_names(Dataset, "../data/pathogen_aliases.csv"),
         Method  = fix_method_names(Method, "../data/method_aliases.csv"),
         Value = ifelse(is.na(Value), 0, Value))

write.table(Xmetrics, 
            "../output/Performance_metrics_all.tsv", sep = "\t",
            quote = FALSE, row.names = FALSE)




## ================ STATISTICAL ANALYSIS ===============
## (1) Internal comparisons: 
## Experimental factor: Method (EpitopeTransfer, NPTransfer, ESM baseline)
## Experimental factor: Embedder (ESM-1b, ESM2)
## Blocking factor: Datasets
## Using blocked factorial ANOVA, Method+Embedder+Dataset
## Note: "Holdout" and "Noleakage" results are the same thing for the internal methods.
##
## (2) External comparisons: 
## EpitopeTransfer(ESM2) vs. external baselines: 
## Experimental factor: Method (EpitopeTransfer(ESM2), Bepipred 3.0, Epidope, EpitopeVec, Epitope1D)
## Blocking factor: Datasets
## Using one-way blocked ANOVA, Method+Dataset
## Main analysis done on No-Leakage sample, but additional analysis also calculated for the internal Holdout one (which has some leakage in terms of performance calculation for external baselines).
##

# Main target metric: 
metr <- "AUC"

# Isolate metrics for the internal tools only
Xmetrint <- Xmetrics %>%
  filter(Set == "NoLeakage", grepl("ESM", Method), Metric == metr) %>%
  rowwise() %>%
  mutate(Embedder = factor(strsplit(as.character(Method), split = "\\(|\\)")[[1]][2]),
         Method   = factor(strsplit(as.character(Method), split = "\\(|\\)")[[1]][1],
                           levels = c("EpitopeTransfer", "NPTransfer", "Baseline"),
                           labels = c("EpitopeTransfer", "NPTransfer", "ESM Baseline"),
                           ordered = TRUE))

# Isolate metrics for the external tools only
Xmetrext <- Xmetrics %>%
  filter(Set == "NoLeakage", !grepl("Baseline|NPT|ESM-1b", Method), Metric == metr) %>%
  mutate(Method = relevel(factor(as.character(Method)), ref = "EpitopeTransfer(ESM2)"))


# Calculate inferential results for internal comparisons
res.internal <- Xmetrics %>% 
  filter(Set == "NoLeakage", grepl("ESM", Method)) %>%
  rowwise() %>%
  mutate(Embedder = factor(strsplit(as.character(Method), split = "\\(|\\)")[[1]][2]),
         Method   = factor(strsplit(as.character(Method), split = "\\(|\\)")[[1]][1],
                           levels = c("EpitopeTransfer", "NPTransfer", "Baseline"),
                           labels = c("EpitopeTransfer", "NPTransfer", "ESM Baseline"),
                           ordered = TRUE)) %>%
  calc_stats2(blockEmbedder = TRUE)

# Calculate inferential results for external comparisons (NoLeakage)
res.external <- Xmetrics %>%
  filter(Set == "NoLeakage", !grepl("Baseline|NPT|ESM-1b", Method)) %>%
  mutate(Method = relevel(factor(as.character(Method)), ref = "EpitopeTransfer(ESM2)")) %>%
  calc_stats2(blockEmbedder = FALSE)

# Estimates for external comparisons (Holdout; for plotting)
res.external2 <- Xmetrics %>%
  filter(Set == "Holdout", !grepl("Baseline|NPT|ESM-1b", Method)) %>%
  mutate(Method = relevel(factor(as.character(Method)), ref = "EpitopeTransfer(ESM2)")) %>%
  calc_stats2(blockEmbedder = FALSE)

# Checking detailed results (internal)
# res.internal[[metr]]$aov.summ
# res.internal[[metr]]$mht.summ
# res.internal[[metr]]$mht.ci
# par(mfrow = c(2,2)); plot(res.internal[[metr]]$aov.model); par(mfrow = c(1,1))

# Checking detailed results (external)
# res.external[[metr]]$aov.summ
# res.external[[metr]]$mht.summ
# res.external[[metr]]$mht.ci
# par(mfrow = c(2,2)); plot(res.external[[metr]]$aov.model); par(mfrow = c(1,1))


# Inference + estimates for full experimental design (internal)
Xmetrint.summ <- Xmetrint %>%
  group_by(Method) %>%
  summarise(Est = mean(Value),
            se  = sd(Value) / sqrt(n()),
            ymin = min(Value),
            .groups = "drop") %>%
  mutate(pValue = c(NA, res.internal[[metr]]$mht.summ$test$pvalues))

# Estimates for results under each embedder (internal; for plotting)
Xmetrint.summ2 <- Xmetrint %>%
  group_by(Method, Embedder) %>%
  summarise(Est = mean(Value),
            se  = sd(Value) / sqrt(n()),
            ymin = min(Value),
            .groups = "drop")


# Inference + estimates for main experiment (external; NoLeakage)
Xmetrext.summ <- Xmetrext %>%
  group_by(Method) %>%
  summarise(Est = mean(Value),
            se  = sd(Value) / sqrt(n()),
            ymin = min(Value),
            .groups = "drop") %>%
  mutate(pValue = c(NA, res.external[[metr]]$mht.summ$test$pvalues))

# Estimates for results (external; Holdout; for plotting)
Xmetrext.summ2 <- Xmetrics %>%
  filter(Set == "Holdout", !grepl("Baseline|NPT|ESM-1b", Method), Metric == metr) %>%
  mutate(Method = relevel(factor(as.character(Method)), ref = "EpitopeTransfer(ESM2)")) %>%
  group_by(Method) %>%
  summarise(Est = mean(Value),
            se  = sd(Value) / sqrt(n()),
            ymin = min(Value),
            .groups = "drop") %>%
  mutate(pValue = c(NA, res.external2[[metr]]$mht.summ$test$pvalues)) %>%
  filter(!grepl("ESM", Method))


## Combined dataset for figure
Xmetrall <- bind_rows(mutate(Xmetrint, Method = as.character(Method)), 
                      mutate(Xmetrext, Method = as.character(Method))) %>%
  mutate(Method = factor(Method,
                         levels = c("EpitopeTransfer", "NPTransfer", 
                                    "ESM Baseline", 
                                    "EpitopeTransfer(ESM2)",
                                    "Bepipred 3.0", "Epidope", 
                                    "Epitope1D", "EpitopeVec"),
                         ordered = TRUE))

## Plot results and inference
Xmetrall%>%
  ggplot(aes(x = Method, y = Value)) + 
  geom_violin(alpha = 0, linewidth = .5, 
              show.legend = FALSE) +
  geom_vline(xintercept = c(3.48, 3.52), colour = "#444444") + 
  geom_jitter(data = Xmetrext, width = 0.1, alpha = .75, colour = "#777777", show.legend = FALSE) + 
  geom_point(data = Xmetrint, 
             aes(colour = Embedder, shape = Embedder),
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.2)) + 
  geom_pointrange(data = Xmetrext.summ2, 
                  aes(x = Method, y = Est, 
                      ymin = Est-se, ymax = Est+se),
                  position = position_nudge(x = 0.05), 
                  size = .75, show.legend = FALSE, shape = 1) +
  geom_pointrange(data = Xmetrext.summ, 
                  aes(x = Method, y = Est, ymin = Est-se, ymax = Est+se),
                  size = 1.25, shape = 15, linewidth = 1,
                  show.legend = FALSE) +
  geom_pointrange(data = Xmetrint.summ, 
                  aes(x = Method, y = Est, ymin = Est-se, ymax = Est+se),
                  size = 1.25, shape = 15, linewidth = 1, 
                  show.legend = FALSE) +
  geom_pointrange(data = Xmetrint.summ2, 
                  aes(x = Method, shape = Embedder, y = Est, 
                      ymin = Est-se, ymax = Est+se, colour = Embedder),
                  position = position_dodge(width = 0.5), size = 1, show.legend = FALSE,
                  alpha = .8,) +
  theme_light() + 
  xlab("") + ylab(metr) +
  # ylim(0.45, 1) +
  annotate("segment",
           x    = c(1, 1, 4, 4, 4, 4),
           xend = c(2, 3, 5, 6, 7, 8),
           y = c(.5, .475, .875, .9, .925, .95),
           yend = c(.5, .475, .875, .9, .925, .95),
           colour = "#777777", linewidth = 1) +
  annotate("label", x = c(4.5, 5, 5.5, 6), y = c(.875, .9, .925, .95), 
           size = 5, label.size = 0,
           label = paste0("italic(P) == ",
                          sprintf("%2.3E", res.external[[metr]]$mht.summ$test$pvalues)),
           parse = TRUE) +
  annotate("label", x = c(1.5, 2), y = c(.5, .475), 
           size = 5, label.size = 0,
           label = paste0("italic(P) == ",
                          sprintf("%2.3E", res.internal[[metr]]$mht.summ$test$pvalues)),
           parse = TRUE) +
  theme(legend.position.inside = c(.25,.9),
        legend.position = "inside",
        legend.background = element_rect(fill = NA, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(hjust = .5, face = "bold", size = 14),
        axis.text.y = element_text(size = 14),
        axis.title  = element_text(face = "bold", size = 14),
        axis.line = element_line(linewidth = .5),
        axis.ticks = element_line(linewidth = .5, colour = "black"),
        #panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) + 
  guides(shape = guide_legend(override.aes = list(size = 3))) + 
  annotate("text", x = 3.3, y = 0.98, label = "(a)", size = 8) +
  annotate("text", x = 8.4, y = 0.98, label = "(b)", size = 8) +
  scale_x_discrete(labels = c("EpitopeTransfer", "NPTransfer", "ESM Baseline",
                              "EpitopeTransfer\n(ESM2)", "Bepipred 3.0", 
                              "Epidope", "Epitope1D", "EpitopeVec")) +
  scale_y_continuous(breaks = .05*(9:20))


ggsave("../figures/all_comparisons.png",
       height = 2000, width = 6000, units = "px")
ggsave("../figures/all_comparisons.pdf",
       height = 2000, width = 6000, units = "px")


# plot gains/losses of EpitopeTransfer(ESM2) vs. baselines for each dataset
gainthres <- 0.05
lossthres <- -0.05
neutralband <- c(-0.05, 0.05)
aliases <- read.csv("../data/method_aliases.csv")

Diff.df <- res.external[[metr]]$df %>%
  group_by(Dataset) %>%
  mutate(Diff = Value[Method == "EpitopeTransfer(ESM2)"] - Value) %>%
  ungroup() %>%
  filter(Method != "EpitopeTransfer(ESM2)") %>%
  mutate(Dataset = ifelse(grepl("\\.", Dataset),
                          paste0("italic(",
                                 gsub(".", "%.%", Dataset, fixed = TRUE),
                                 ")"),
                          Dataset))

Diff.summ <- Diff.df  %>%
  group_by(Method) %>%
  summarise(Ngain = sum(Diff > 0),
            Nloss = sum(Diff < 0),
            GainRatio = 100 * Ngain / n(),
            meanGain = mean(Diff[Diff > 0]),
            meanLoss = mean(Diff[Diff < 0]),
            seGain   = sd(Diff[Diff > 0]) / sqrt(Ngain),
            seLoss   = sd(Diff[Diff < 0]) / sqrt(Nloss),
            .groups = "drop") %>%
  mutate(across(starts_with("se"), ~ifelse(is.na(.x), 0, .x))) %>%
  ungroup() %>% 
  pivot_longer(meanGain:meanLoss, names_to = "type", values_to = "mean") %>%
  pivot_longer(seGain:seLoss, names_to = "setype", values_to = "se") %>%
  pivot_longer(Ngain:Nloss, names_to = "gain", values_to = "N") %>%
  mutate(type = gsub("mean", "", type),
         setype = gsub("se", "", setype),
         gain   = ifelse(gain == "Ngain", "Gain", "Loss")) %>%
  filter(type == setype,
         type == gain) %>%
  dplyr::select(Method, type, mean, se, N)

Diff.summ %>%
  ggplot(aes(x = Method, y = mean)) + 
  geom_rect(xmin = 0, xmax = 5.5, ymin = neutralband[1], ymax = neutralband[2],
            colour = NA, fill = "#cccccc22") +
  geom_violinhalf(data = Diff.df %>% filter(Diff > 0), 
                  aes(x = Method, y = Diff, group = Method), 
                  flip = TRUE, scale = "count", linewidth = NA,
                  fill = "#bbffbbaa") +
  geom_violinhalf(data = Diff.df %>% filter(Diff < 0), 
                  aes(x = Method, y = Diff, group = Method), 
                  flip = FALSE, scale = "count", linewidth = NA,
                  fill = "#ffbbbbaa") +
  geom_pointrange(aes(ymin = mean - se, ymax = mean + se, group = type), 
                  position = position_dodge(width = .15), 
                  size = 1, linewidth = 1, 
                  colour = "#44444455", stroke = 0,
                  pch = 15) +
  geom_text_repel(data = Diff.df %>% group_by(Dataset) %>% filter(all(Diff > gainthres)),
                  aes(x = Method, y = abs(Diff), label = Dataset),
                  nudge_x = .1,
                  direction = "y",
                  hjust = 0,
                  min.segment.length = 0, force = 5,
                  segment.colour = "#66666699",
                  size = 4.5, box.padding = 0.5,
                  show.legend = FALSE, parse = TRUE) +
  geom_text_repel(data = Diff.df %>% filter(Diff < min(neutralband)),
                  aes(x = Method, y = Diff, label = Dataset),
                  size = 4.5, , box.padding = 0.5,
                  min.segment.length = 0,  force = 5,
                  nudge_x = -0.1,
                  direction = "y",
                  hjust = 1,
                  show.legend = FALSE,
                  segment.colour = "#66666699",
                  colour = "#888888", parse = TRUE) +
  geom_point(data = Diff.df %>% filter(Diff > 0), 
             aes(x = Method, y = abs(Diff), colour = (Diff > gainthres)),
             size = 2, pch = 20,
             show.legend = FALSE) + 
  geom_point(data = Diff.df %>% filter(Diff < 0),
             aes(x = Method, y = Diff, colour = (Diff < -0.2)),
             size = 2, pch = 20,
             show.legend = FALSE) +
  theme_light() + 
  theme(axis.text.x = element_text(hjust = .5, face = "bold", size = 14),
        axis.text.y = element_text(size = 14),
        axis.title  = element_text(face = "bold", size = 14),
        axis.line = element_line(linewidth = .5),
        axis.ticks = element_line(linewidth = .5, colour = "black")) + 
  xlab("") + ylab("EpitopeTransfer gains in AUC") +
  scale_colour_manual(values = c("#777777", "#ff0000")) + 
  scale_x_discrete(labels = c("EpitopeTransfer(ESM2)\nvs. Bepipred 3.0", 
                              "EpitopeTransfer(ESM2)\nvs. Epidope", 
                              "EpitopeTransfer(ESM2)\nvs. Epitope1D", 
                              "EpitopeTransfer(ESM2)\nvs. EpitopeVec")) +
  geom_hline(yintercept = 0) +
  annotate("text", x = 4.5, y = 0.425, label = "(c)", size = 8) +
  scale_y_continuous(limits = c(-0.15, 0.45), breaks = 0.05*(-3:9), expand = c(0.001,0))

ggsave(paste0("../figures/winloss_violin_AUC.png"), 
       width = 6000, height = 2000, units = "px")
ggsave(paste0("../figures/winloss_violin_AUC.pdf"), 
       width = 6000, height = 2000, units = "px")


# Save relevant results in a more human readable format
stats.table.int <- bind_rows(lapply(res.internal, \(x) x$summary))
stats.table.ext <- bind_rows(lapply(res.external, \(x) x$summary))

rownames(stats.table.int) <- NULL
rownames(stats.table.ext) <- NULL

write.table(stats.table.int, 
            "../output/comparisons_internal.tsv", sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(stats.table.ext, 
            "../output/comparisons_external.tsv", sep = "\t",
            quote = FALSE, row.names = FALSE)

write.table(Xmetrics, "../output/all_results.tsv", sep = "\t",
            quote = FALSE, row.names = FALSE)


## Extract some further info for paper discussion
Diff.df %>%
  group_by(Dataset) %>%
  mutate(Sel = all(Diff > 0)) %>%
  filter(Sel) %>%
  summarise(Sel = first(Sel))

Diff.df %>%
  group_by(Dataset) %>%
  mutate(Sel = all(Diff > 0.05)) %>%
  filter(Sel) %>%
  summarise(Sel = first(Sel))

Diff.df %>%
  group_by(Dataset) %>%
  mutate(Sel = all(Diff < 0)) %>%
  filter(Sel) %>%
  arrange(desc(Diff)) %>%
  print(n = 100)


Diff.df %>%
  filter(Dataset %in% c("S. mansoni", "M. tuberculosis", "O. volvulus", "B. pertussis", "Orthopoxvirus"))

