library(dplyr)
library(tidyr)
library(parallel)
library(pbapply)
source("myfuns.R")

ncpus <- 10
cl = parallel::makePSOCKcluster(ncpus)

# # Process Epitope1D results (only needs to be done once)
# cat("\nProcessing Epitope1D results\n")
# process_epitope1d(resdir = "../results/epitope1d/",
#                   datadir = "../data/lindeberg/",
#                   cl = cl)

# Load and consolidate results
fp <- dir("../results", pattern = "prob\\_by", full.names = TRUE, recursive = TRUE)

cat("\nReading results\n")
X1 <- pblapply(fp[!grepl("esm|ESM", fp)],
              function(myfile){
                parts <- strsplit(myfile, split = "/")[[1]]
                x <- read.csv(myfile, header = TRUE)
                x$Embedder  <- "other"
                x$Predictor <- gsub("\ ", "\\_", tolower(parts[3]))
                x$Pathogen  <- gsub("\ ", "\\_", tolower(parts[4]))
                x$Fold      <- parts[5]
                return(x)
              }, cl = cl) %>%
  bind_rows() %>%
  rename(Info_protein_id = Protein_id,
         Info_pos        = Position,
         Class           = True_class,
         Prob            = Predicted_prob)

## To check for consistency
# table(X1$Pathogen, X1$Predictor)

X2 <- pblapply(fp[grepl("esm|ESM", fp)],
               function(myfile){
                 parts <- strsplit(myfile, split = "/")[[1]]
                 x <- read.csv(myfile, header = TRUE)
                 x$Embedder  <- tolower(gsub("-", "", parts[3]))
                 x$Predictor <- gsub("\ ", "\\_", tolower(parts[4]))
                 x$Pathogen  <- gsub("\ ", "\\_", tolower(parts[5]))
                 x$Fold      <- parts[6]
                 return(x)
               }, cl = cl) %>%
  bind_rows() %>%
  rename(Info_protein_id = Protein_id,
         Info_pos        = Position,
         Class           = True_class,
         Prob            = Predicted_prob)

## To check for consistency
# table(X2$Pathogen, X2$Predictor)

X <- bind_rows(X1, X2)

## To check for consistency
# table(X$Pathogen, X$Predictor, X$Embedder)

# Load and prepare data about the external predictors' training sets
inTraining <- readRDS("../data/peptides_in_baseline_training_sets.rds") %>%
  rowwise() %>%
  mutate(ll = Info_end_pos - Info_start_pos + 1) %>%
  uncount(ll) %>%
  group_by(Info_PepID) %>%
  mutate(Info_pos = Info_start_pos + 1:n() - 1) %>%
  ungroup() %>%
  dplyr::select(-Info_start_pos, -Info_end_pos, -Info_sequence, -Info_peptide, -InTraining_BP2, -Info_PepID) %>%
  dplyr::select(Info_protein_id, Info_pos, everything()) %>%
  pivot_longer(InTraining_BP3:InTraining_Ep1D, names_to = "Predictor", values_to = "InTraining")

inTraining$Predictor[grep("BP3", inTraining$Predictor)] <- "bepipred3"
inTraining$Predictor[grep("EpDop", inTraining$Predictor)] <- "epidope"
inTraining$Predictor[grep("EpVec", inTraining$Predictor)] <- "epitopevec"
inTraining$Predictor[grep("Ep1D", inTraining$Predictor)] <- "epitope1d"

inTraining <- inTraining %>% filter(InTraining == 1, 
                                    Info_protein_id %in% X$Info_protein_id)

X <- X %>%
  left_join(inTraining, by = c("Info_protein_id", "Info_pos", "Predictor")) %>%
  mutate(InTraining = ifelse(is.na(InTraining), 0, InTraining))

# Incorporate fold information into Epitope1D results
a <- X %>% filter(Predictor == "bepipred3") %>% dplyr::select(Pathogen, Info_protein_id, Info_pos, Fold)
b <- X %>% filter(Predictor == "epitope1d") %>% dplyr::select(-Fold) %>%
  left_join(a, by = c("Pathogen", "Info_protein_id", "Info_pos"))
Y <- X %>%
  filter(Predictor != "epitope1d") %>%
  bind_rows(b)


saveRDS(Y, "../data/results/consolidated_results_v9.rds")

parallel::stopCluster(cl)
gc()
