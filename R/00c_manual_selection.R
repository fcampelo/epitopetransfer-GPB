library(dplyr)
library(epitopes)

epitopes <- readRDS("../data/epitopes.rds")
proteins <- readRDS("../data/proteins.rds")
taxonomy <- readRDS("../data/taxonomy.rds")

ncpus <- max(1, min(20, parallel::detectCores() - 1))

tmp <- readRDS("../data/orgdata/txIDs_filtered.rds")
txIDs <- tmp$txIDs
myOrgs <- tmp$myOrgs
rm(tmp)
txIDs <- txIDs[-22]
myOrgs <- myOrgs[-22, ]

for (i in seq_along(txIDs)){
  cat("\ni = ", i)
  fn <- paste0("../data/orgdata/",
               txIDs[[i]]$UID[1], "_for_", txIDs[[i]]$UID[nrow(txIDs[[i]])],
               "_", txIDs[[i]]$ScientificName[nrow(txIDs[[i]])], ".rds")
  x <- readRDS(fn)
  tmp <- x$peptides %>%
    group_by(Info_protein_id) %>%
    summarise(Info_organism_id = first(Info_organism_id),
              Info_cluster = first(Info_cluster),
              Info_split   = first(Info_split))
  
  prots <- x$proteins
  prots <- prots  %>%
    left_join(tmp, by = c("UID" = "Info_protein_id")) %>%
    select(Info_organism_id,
           Info_protein_id = UID,
           Info_sequence = TSeq_sequence,
           Info_cluster,
           Info_split)
  
  txIDs[[i]]$IDs <- ""
  txIDs[[i]]$IDs[1] <- paste(unique(x$df$Info_organism_id), collapse = ",")
  cat(":")
  for (j in 2:nrow(txIDs[[i]])){
    tmp <- x$df %>%
      filter_epitopes(tax_list = taxonomy, orgIDs = txIDs[[i]]$UID[j], 
                      orgID_column = "Info_organism_id", 
                      hostID_column = "Info_host_id")
    txIDs[[i]]$IDs[j] <- paste(unique(tmp$Info_organism_id), collapse = ",")
    cat(".")
  }
  
  dn <- paste0("../data/lindeberg/", 
               gsub("\\.\\ ", "", myOrgs$Lowest.Name[i]))
  if(!dir.exists(dn)) dir.create(dn, recursive = TRUE)
  saveRDS(prots,
          paste0(dn, "/proteins_", txIDs[[i]]$UID[1], "_for_", 
                 txIDs[[i]]$UID[nrow(txIDs[[i]])], "_", 
                 txIDs[[i]]$ScientificName[nrow(txIDs[[i]])], ".rds"))
  saveRDS(x$peptides,
          paste0(dn, "/peptides_", txIDs[[i]]$UID[1], "_for_", 
                 txIDs[[i]]$UID[nrow(txIDs[[i]])], "_", 
                 txIDs[[i]]$ScientificName[nrow(txIDs[[i]])], ".rds"))
  saveRDS(txIDs[[i]], 
          paste0(dn, "/taxonomy_", txIDs[[i]]$UID[1], "_for_", 
                 txIDs[[i]]$UID[nrow(txIDs[[i]])], "_", 
                 txIDs[[i]]$ScientificName[nrow(txIDs[[i]])], ".rds"))
  
  write.csv(prots,
          paste0(dn, "/proteins_", txIDs[[i]]$UID[1], "_for_", 
                 txIDs[[i]]$UID[nrow(txIDs[[i]])], "_", 
                 txIDs[[i]]$ScientificName[nrow(txIDs[[i]])], ".csv"), 
          row.names = FALSE)
  write.csv(x$peptides,
          paste0(dn, "/peptides_", txIDs[[i]]$UID[1], "_for_", 
                 txIDs[[i]]$UID[nrow(txIDs[[i]])], "_", 
                 txIDs[[i]]$ScientificName[nrow(txIDs[[i]])], ".csv"), 
          row.names = FALSE)
  write.csv(txIDs[[i]], 
          paste0(dn, "/taxonomy_", txIDs[[i]]$UID[1], "_for_", 
                 txIDs[[i]]$UID[nrow(txIDs[[i]])], "_", 
                 txIDs[[i]]$ScientificName[nrow(txIDs[[i]])], ".csv"), 
          row.names = FALSE)
  
}


