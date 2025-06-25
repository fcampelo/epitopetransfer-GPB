library(dplyr)
library(epitopes)

epitopes <- readRDS("../data/epitopes.rds")
proteins <- readRDS("../data/proteins.rds")
taxonomy <- readRDS("../data/taxonomy.rds")

ncpus <- max(1, min(20, parallel::detectCores() - 1))


# Organism IDs:
myOrgs <- matrix(c("HIV-1", "Virus", "11676",
                   "SARS-CoV-2", "Virus", "694009",
                   "N. gonorrhoeae", "Bacteria", "485",
                   "M. tuberculosis", "Bacteria", "1773",
                   "P. aeruginosa", "Bacteria", "287",
                   "S. enterica", "Bacteria", "28901",
                   "B. pertussis", "Bacteria", "520",
                   "C. diphtheriae", "Bacteria", "1717",
                   "E. coli", "Bacteria", "562",
                   "C. neoformans", "Fungi", "5207",
                   "P. falciparum", "Protozoa", "5833",
                   "G. intestinalis", "Protozoa", "5741",
                   "T. gondii", "Protozoa", "5811",
                   "A. lumbricoides", "Helminth", "6252",
                   "S. mansoni", "Helminth", "6183"),
                 ncol = 3, byrow = TRUE) %>%
  as.data.frame() %>%
  rename(Name = V1, Type = V2, ID = V3)

txIDs <- get_taxonomy(myOrgs$ID)
txIDs <- mapply(
  function(x, y) {
    a <- rbind(x$Taxonomy, c(y, "organism", x$UID))
    a$nPos <- 0
    a$nNeg <- 0
    a$nTot <- 0
    a[a$UID != "131567", ]
  },
  x = txIDs, y = myOrgs$Name,
  SIMPLIFY = FALSE)

for (i in seq_along(txIDs)){
  for (j in seq_along(txIDs[[i]]$UID)){
    cat("\nProcessing ", i, ":", j, "/", length(txIDs[[i]]$UID), "\n")
    fn <- paste0("../data/orgdata/", txIDs[[i]]$UID[j], ".rds")
    if(j == 1) {
      if(file.exists(fn)){
        topData <- readRDS(fn)
      } else {
        topData <- epitopes %>%
          filter_epitopes(orgIDs        = txIDs[[i]]$UID[j],
                          orgID_column  = "sourceOrg_id",
                          hostID_column = "host_id",
                          tax_list      = taxonomy) %>%
          consolidate_data(proteins, only_exact = FALSE, ncpus = ncpus) %>%
          extract_peptides(min_peptide = 5, max_epitope = 30, window_size = 15)
        
        topData$proteins <- proteins %>%
          dplyr::filter(UID %in% unique(topData$peptides$Info_protein_id)) %>%
          dplyr::select(TSeq_taxid, UID, TSeq_sequence) %>%
          dplyr::rename(Info_organism_id = TSeq_taxid,
                        Info_sequence = TSeq_sequence)
        
        saveRDS(topData, fn)
      }
      txIDs[[i]]$nTot[j] <- nrow(topData$peptides)
      txIDs[[i]]$nPos[j] <- sum(topData$peptides$Class == 1)
      txIDs[[i]]$nNeg[j] <- sum(topData$peptides$Class == -1)
      
    } else {
      if(file.exists(fn)){
        lowData <- readRDS(fn)
      } else {
        tmp <- epitopes %>%
          filter_epitopes(orgIDs        = txIDs[[i]]$UID[j],
                          orgID_column  = "sourceOrg_id",
                          hostID_column = "host_id",
                          tax_list      = taxonomy)
        tmp <- unique(tmp$sourceOrg_id)
        lowData <- topData
        lowData$df <- dplyr::filter(lowData$df,
                                    Info_organism_id %in% tmp)
        lowData$peptides <- dplyr::filter(lowData$peptides,
                                          Info_organism_id %in% tmp)
        saveRDS(lowData, fn)
      }
      txIDs[[i]]$nTot[j] <- nrow(lowData$peptides)
      txIDs[[i]]$nPos[j] <- sum(lowData$peptides$Class == 1)
      txIDs[[i]]$nNeg[j] <- sum(lowData$peptides$Class == -1)
    }
  }
}
rm(topData, lowData)

saveRDS(txIDs, "../data/orgdata/txIDs.rds")


# Remove datasets with too little low-level data
txIDs[c(3, 10, 12, 14)] <- NULL
myOrgs <- myOrgs[-c(3, 10, 12, 14), ]


# Define lowest level for each case
txIDs[[1]] <- txIDs[[1]][1:9, ]
txIDs[[5]] <- txIDs[[5]][1:5, ]
txIDs[[7]] <- txIDs[[7]][1:7, ]
myOrgs$Lowest.Name <- sapply(txIDs, function(x) x$ScientificName[nrow(x)])
myOrgs$Lowest.ID   <- sapply(txIDs, function(x) x$UID[nrow(x)])

saveRDS(list(txIDs  = txIDs, 
             myOrgs = myOrgs),
        "../data/orgdata/txIDs_filtered.rds")


# Define relevant data splits
diss_matrix <- readRDS("../data/protein_dissimilarity.rds")
for (i in seq_along(txIDs)){
  cat("\n\n\n\n\n\ni = ", i, "\n\n")
  x <- readRDS(paste0("../data/orgdata/", txIDs[[i]]$UID[1], ".rds"))
  ids <- x$df %>%
    filter_epitopes(tax_list = taxonomy,
                    orgIDs = myOrgs$Lowest.ID[i],
                    orgID_column = "Info_organism_id",
                    hostID_column = "Info_host_id")
  ids <- unique(ids$Info_organism_id)
  x <- x %>%
    epitopes::make_data_splits(proteins,
                               split_level = "protein",
                               target_props = rep(.2, 5),
                               similarity_threshold = .7,
                               diss_matrix = diss_matrix,
                               id_force_splitting = ids,
                               tax_list = taxonomy,
                               ncpus = ncpus)
  saveRDS(x, paste0("../data/orgdata/",
                    txIDs[[i]]$UID[1], "_for_", txIDs[[i]]$UID[nrow(txIDs[[i]])],
                    "_", txIDs[[i]]$ScientificName[nrow(txIDs[[i]])], ".rds"))
}


