# This script requires DIAMOND to be installed.
# Please download the appropriate version for your OS from
# https://github.com/bbuchfink/diamond/releases/tag/v2.0.15
#
# and place the un-compressed files in folder code/diamond/

library(dplyr)
library(tidyr)
library(seqinr)

proteins <- readRDS("../data/proteins.rds")

seqinr::write.fasta(sequences = as.list(proteins$TSeq_sequence), 
                    names     = proteins$UID, 
                    file.out  = "../data/proteins.fa")

myseqs <- proteins$TSeq_sequence
myseqs <- gsub("B|J|O|U|X|Z", "<mask>", myseqs)
seqinr::write.fasta(sequences = as.list(myseqs), 
                    names     = proteins$UID, 
                    file.out  = "../data/proteins_masked.fa")

# Get dissimilarity matrix and hierarchical clustering of the proteins:
# Run DIAMOND to get local alignment scores
message("Calculating similarities using DIAMOND...")

# Run DIAMOND
# create a diamond-formatted database file
if(!dir.exists("../data/diamond")) dir.create("../data/diamond", recursive = TRUE)
system("diamond makedb --in ../data/proteins.fa -d ../data/diamond/proteins-reference")

# running a search in blastp mode
system("diamond blastp -d ../data/diamond/proteins-reference -q ../data/proteins.fa --ultra-sensitive -o ../data/protein-matches.tsv")

toRM <- dir("./", pattern = ".tmp")
if(length(toRM) > 0) ignore <- file.remove(toRM)

# Read results
scores <- read.csv("../data/protein-matches.tsv", sep = "\t",
                  header = FALSE,
                  stringsAsFactors = FALSE)
names(scores) <- c("qseqid", "sseqid", "pident", "length",
                   "mismatch", "gapopen", "qstart", "qend",
                   "sstart", "send", "evalue", "bitscore")

scores <- scores %>%
  filter(length >= 8) %>%
  group_by(qseqid, sseqid) %>%
  arrange(desc(pident)) %>%
  summarise(across(everything(), first), .groups = "drop") %>%
  mutate(diss = 1 - pident/100) %>%
  select(qseqid, sseqid, diss) %>%
  pivot_wider(names_from = sseqid, values_from = diss, values_fill = NA) %>%
  arrange(qseqid) %>%
  as.data.frame()

rownames(scores) <- scores$qseqid

# Make the dissimilarity scores matrix
scores <- scores %>%
  select(order(colnames(scores)),
         -qseqid)

pnames <- rownames(scores)

protIDs <- names(seqinr::read.fasta("../data/proteins.fa", as.string = TRUE))
missing <- protIDs[which(!(protIDs %in% rownames(scores)))]

scores[(nrow(scores) + 1):(nrow(scores) + length(missing)), ] <- 1
scores[, (ncol(scores) + 1):(ncol(scores) + length(missing))] <- 1

scores <- as.matrix(scores)
diag(scores) <- 0

colnames(scores) <- c(pnames, missing)
rownames(scores) <- c(pnames, missing)
scores[is.na(scores)] <- 1

saveRDS(scores, "../data/protein_dissimilarity.rds")
