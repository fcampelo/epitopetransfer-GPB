library(dplyr)
library(seqinr)

dl <- dir("../data/lindeberg/", full.names = TRUE)

for (i in seq_along(dl)){
  proteins <- readRDS(dir(dl[i], pattern = "^protein.+\\.rds", full.names = TRUE))
  peptides <- readRDS(dir(dl[i], pattern = "^pept.+\\.rds", full.names = TRUE))
  taxonomy <- readRDS(dir(dl[i], pattern = "^taxon.+\\.rds", full.names = TRUE))
  
  proteins <- proteins %>%
    filter(Info_organism_id %in% strsplit(taxonomy$IDs[nrow(taxonomy)], split = ",")[[1]])
  
  peptides <- peptides %>%
    filter(Info_organism_id %in% strsplit(taxonomy$IDs[nrow(taxonomy)], split = ",")[[1]])
  
  write.fasta(sequences = as.list(proteins$Info_sequence), 
              names = proteins$Info_protein_id, 
              file.out = paste0(dl[i], "/LowerLevelProts.fa"))
  
  write.fasta(sequences = as.list(peptides$Info_peptide), 
              names = peptides$Info_PepID, 
              file.out = paste0(dl[i], "/LowerLevelPepts.fa"))
}