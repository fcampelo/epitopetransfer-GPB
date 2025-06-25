library(seqinr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(boot)
library(tidyr)
library(rlang)

fl <- dir("../data/lindeberg/", full.names = TRUE)

for (i in seq_along(fl)){
  tmp <- readRDS(dir(fl[i], full.names = TRUE)[grep("peptides.+\\.rds$", dir(fl[i]))]) %>%
    dplyr::left_join(readRDS(dir(fl[i], full.names = TRUE)[grep("protein.+\\.rds", dir(fl[i]))]), 
                     by = "Info_protein_id") %>%
    dplyr::select(Info_PepID, Info_protein_id, Info_start_pos, Info_end_pos, Info_peptide, Info_sequence)
  if (i == 1) df <- tmp else df <- rbind(df, tmp)
}

df <- df %>%
  group_by(Info_peptide, Info_protein_id) %>%
  summarise(across(everything(), first))

prots <- df %>%
  dplyr::group_by(Info_protein_id) %>%
  dplyr::summarise(Info_sequence = dplyr::first(Info_sequence))

seqinr::write.fasta(sequences = as.list(df$Info_peptide),
                    names     = df$Info_PepID,
                    file.out  = "../data/unique_peptides.fa")
seqinr::write.fasta(sequences = as.list(prots$Info_sequence),
                    names     = prots$Info_protein_id,
                    file.out  = paste0("../data/unique_prots.fa"))

# ==============================
# Flag which  entries are included in the training sets of other predictors

thres <- 0.8
preds <- c("BP2", "BP3", "EpDop", "EpVec", "Ep1D")

for (i in seq_along(preds)){
  system(paste0("diamond/diamond makedb --in ../data/baselines/", preds[i], "_training.fa -d ../data/diamond/", preds[i], "_entries"))
  
  
  outfile <- paste0("../data/baselines/", preds[i], "-matches.tsv")
  system(paste0("diamond/diamond blastp -d ../data/diamond/", preds[i], "_entries ", 
                "-q ../data/unique_peptides.fa --ultra-sensitive -o ", outfile))

  scores <- read.table(outfile, sep = "\t",
                       header = FALSE, stringsAsFactors = FALSE)
  names(scores) <- c("qseqid", "sseqid", "pident", "length",
                     "mismatch", "gapopen", "qstart", "qend",
                     "sstart", "send", "evalue", "bitscore")
  
  scores <- scores %>% 
    filter(pident > 100*thres) %>% 
    rename(Info_PepID = qseqid) %>%
    group_by(Info_PepID) %>%
    arrange(desc(pident)) %>%
    summarise(across(everything(), ~ first(.x))) %>%
    mutate(XX = 1) %>%
    dplyr::select(Info_PepID, XX)
  names(scores)[2] <- paste0("InTraining_", preds[i])
  
  df <- df %>%
    left_join(scores, by = "Info_PepID")
}

df <- df %>%
  mutate(across(everything(), ~ifelse(is.na(.x), 0, .x)))

saveRDS(df, "../data/peptides_in_baseline_training_sets.rds")
#===============================================================================
