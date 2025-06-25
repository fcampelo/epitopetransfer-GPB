library(dplyr)
library(epitopes)

epitopes <- readRDS("../data/epitopes.rds")
proteins <- readRDS("../data/proteins.rds")
taxonomy <- readRDS("../data/taxonomy.rds")

ncpus <- max(1, min(20, parallel::detectCores() - 4))


# Organism IDs:
myOrgs <- matrix(c(
  #   "HIV-1", "Virus", "11676",
  #   "SARS-CoV-2", "Virus", "694009",
  #   "N. gonorrhoeae", "Bacteria", "485",
  #   "M. tuberculosis", "Bacteria", "1773",
  #   "P. aeruginosa", "Bacteria", "287",
  #   "S. enterica", "Bacteria", "28901",
  #   "B. pertussis", "Bacteria", "520",
  #   "C. diphtheriae", "Bacteria", "1717",
  #   "E. coli", "Bacteria", "562",
  #   "C. neoformans", "Fungi", "5207",
  #   "P. falciparum", "Protozoa", "5833",
  #   "G. intestinalis", "Protozoa", "5741",
  #   "T. gondii", "Protozoa", "5811",
  #   "A. lumbricoides", "Helminth", "6252",
  #   "S. mansoni", "Helminth", "6183",
    "Orthopoxvirus", "Virus", "10242",
    "Influenza A", "Virus", "11320",
    "Measles morbilivirus", "Virus", "11234",
    "Dengue virus", "Virus", "12637",
    "O. volvulus", "Nemathode", "6282",
    "H. hominis", "Virus", "3052230",
    "Human Gammaherpesvirus 4", "Virus", "10376",
    "R. rickettsii", "Bacteria", "783",
    "S. aureus", "Bacteria", "1280",
    "C. jejuni", "Bacteria", "197",
    "C. trachomatis", "Bacteria", "813",
    "Leptospira", "Bacteria", "171",
    "S. pyogenes", "Bacteria", "1314",
    "K. pneumoniae", "Bacteria", "573",
    "A. baumannii", "Bacteria", "470",
    "Y. pestis", "Bacteria", "632",
    "H. nipahense", "Virus", "3052225",
    "T. brucei", "Protozoa", "5691",
    "L. infantum", "Protozoa", "5671",
    "O. zairense", "Virus", "3052462",
    "P. mirabilis", "Bacteria", "584",
    "H. influenzae", "Bacteria", "727",
    "S. marcescens", "Bacteria", "615",
    "M. pneumoniae", "Bacteria", "2104",
    "L. pneumophila", "Bacteria", "446",
    "C. difficile", "Bacteria", "1496",
    "L. monocytogenes", "Bacteria", "1639",
    "L. rabies", "Virus", "11292",
    "Orthohantavirus", "Virus", "1980442",
    "Alphavirus", "Virus", "11019",
    "Yellow fever virus", "Virus", "11089",
    "L. mammarenavirus", "Virus", "3052310",
    "Enterovirus C", "Virus", "138950 "),
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
    
    rmIDs <- c("131567")   # Supergroup Cellular organisms
    a <- a[!(a$UID %in% rmIDs), ] # Remove supergroups, superkingdoms etc. above what would be useful for this experiment.
    return(a)
  },
  x = txIDs, y = myOrgs$Name,
  SIMPLIFY = FALSE)

for (i in seq_along(txIDs)){
  #cat("\n***** Processing", i, "of", length(txIDs), "*****")
  for (j in seq_along(txIDs[[i]]$UID)){
    cat("\nProcessing ", i, ":", j, "/", length(txIDs[[i]]$UID), "\n")
    fn <- paste0("../data/orgdata/", txIDs[[i]]$UID[j], ".rds")
    if(j == 1) {
      if(file.exists(fn)){
        #    cat("\nReading from file: ",  txIDs[[i]]$ScientificName[j])
        topData <- readRDS(fn)
      } else {
        #    cat("\nGenerating file for: ",  txIDs[[i]]$ScientificName[j])
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
        #   cat("\n\tReading from file: ",  txIDs[[i]]$ScientificName[j])
        lowData <- readRDS(fn)
      } else {
        #  cat("\n\tGenerating file for: ",  txIDs[[i]]$ScientificName[j])
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
# Selection criteria: 
# * at least 10 peptides in each class, at least 20 peptides total in one of the lowest levels (not necessarily the lowest one).
# * At least 50 peptides in each class, at least 150 peptides total in one of the higher levels (below Superkingdom) from taxa distinct from the lower-level one
toRM <- c(10, 12, 15, 18, 21, 24, 25, 27, 32)


txIDs <- txIDs[-toRM]
myOrgs <- myOrgs[-toRM, ]


# Update lowest and highest level where needed and remove empty rows
txIDs <- lapply(txIDs, function(x) x[(x$nTot < 4000) & (x$nTot > 0), ])
txIDs[[1]]  <- txIDs[[1]][-1, ]
txIDs[[8]]  <- txIDs[[8]][1:5, ]
txIDs[[10]]  <- txIDs[[10]][7:8, ]
txIDs[[11]]  <- txIDs[[11]][-1, ]
txIDs[[14]]  <- txIDs[[14]][2:6, ]
txIDs[[15]]  <- txIDs[[15]][6:11, ]
txIDs[[19]]  <- txIDs[[19]][-1, ]
txIDs[[20]]  <- txIDs[[20]][1:4, ]

myOrgs$Lowest.Name <- sapply(txIDs, function(x) x$ScientificName[nrow(x)])
myOrgs$Lowest.ID   <- sapply(txIDs, function(x) x$UID[nrow(x)])


saveRDS(list(txIDs  = txIDs, 
             myOrgs = myOrgs),
        "../data/orgdata/txIDs_filtered.rds")


# Define relevant data splits
diss_matrix <- readRDS("../data/protein_dissimilarity.rds")
for (i in seq_along(txIDs)){
  cat("\nRunning", i, "of", length(txIDs))
  fn <-  paste0("../data/orgdata/",
                txIDs[[i]]$UID[1], "_for_", txIDs[[i]]$UID[nrow(txIDs[[i]])],
                "_", txIDs[[i]]$ScientificName[nrow(txIDs[[i]])], ".rds")
  
  if(!file.exists(fn)){
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
                                 SAopts = list(torun = FALSE),
                                 ncpus = ncpus)
    saveRDS(x, fn)
  }
}


