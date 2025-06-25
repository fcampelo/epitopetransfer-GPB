process_epitope1d <- function(resdir, datadir, cl = 1){
  ### Process Epitope1D predictions
  ### CSV files with predictions was retrieved from Epitope1D in September 2024.
  
  dirs         <- dir(resdir, full.names = TRUE)
  dirnames     <- dir(resdir, full.names = FALSE)
  datadirs     <- dir(datadir, full.names = TRUE)
  datadirnames <- dir(datadir, full.names = FALSE)
  
  mymean <- function(z){ifelse(all(is.na(z)), NA, mean(z, na.rm = TRUE))}
  mymax <- function(z){ifelse(all(is.na(z)), NA, max(z, na.rm = TRUE))}
  
  for (i in seq_along(dirs)){
    preds <- read.csv(dir(dirs[i], pattern = "output\\_prediction.csv", full.names = TRUE), header = TRUE) %>%
      dplyr::select(Info_protein_id = Fasta_header, 
                    Info_Peptide = Peptide, Prob = Score_Epitope)
    
    df    <- readRDS(dir(datadirs[i], pattern = "^peptides.+\\.rds", full.names = TRUE)) %>%
      dplyr::select(Info_PepID, Info_protein_id, Info_peptide, Info_start_pos, Info_end_pos, Class)
    
    prots <- readRDS(dir(datadirs[i], pattern = "^proteins.+\\.rds", full.names = TRUE)) %>%
      dplyr::select(Info_protein_id, Info_sequence) %>%
      dplyr::filter(Info_protein_id %in% unique(preds$Info_protein_id))
    
    cat(sprintf("\nProcessing dir %03d/%03d (%s | %s): %03d proteins\n",
                i, length(dirs),
                dirnames[i], datadirnames[i],
                nrow(prots)))
    
    # Aggregate Epitope1D data by residue.
    x <- pbapply::pblapply(seq_along(prots$Info_protein_id), 
                           function(j, prots, df, preds){
                             require(dplyr)
                             myseq <- prots$Info_sequence[j]
                             tmp   <- preds %>%
                               filter(Info_protein_id == prots$Info_protein_id[j])
                             tmp <- tmp %>%
                               bind_cols(stringr::str_locate(myseq, tmp$Info_Peptide)) %>%
                               filter(!is.na(start))
                             
                             labels <- df %>% 
                               filter(Info_protein_id == prots$Info_protein_id[j]) %>%
                               mutate(N = nchar(Info_peptide)) %>%
                               tidyr::uncount(N) %>%
                               group_by(Info_PepID) %>%
                               mutate(Position = Info_start_pos + (1:n()) - 1) %>%
                               ungroup() %>%
                               group_by(Position) %>%
                               summarise(True_class = sign(sum(Class)), .groups = "drop") %>%
                               mutate(True_class = ifelse(True_class == 0, NA, True_class))
                             
                             x <- data.frame(Pathogen  = "", 
                                             Protein_id = prots$Info_protein_id[j],
                                             Position = 1:nchar(myseq),
                                             Predictor = "", 
                                             Predicted_prob = NA) %>%
                               left_join(labels, by = "Position") %>%
                               rowwise() %>%
                               mutate(Predicted_prob = mymax(tmp$Prob[tmp$start <= Position & tmp$end >= Position]))
                             
                             return(x)
                           }, 
                           prots = prots, df = df, preds = preds, cl = cl) %>%
      bind_rows() %>%
      dplyr::select(Pathogen, Protein_id, Position,
                    True_class, Predictor, Predicted_prob)
    
    write.csv(x, paste0(dirs[i], "/prob_by_position.csv"), 
              row.names = FALSE, quote = FALSE)
    
    
  }
  
  invisible(TRUE)
}



check_seqtype <- function(seq){
  unique_chars <- unique(strsplit(toupper(seq), split = "")[[1]])
  
  idx <- which(!(unique_chars %in% LETTERS))
  if(length(idx) > 0) unique_chars <- unique_chars[-idx]
  
  if(all(unique_chars %in% c("A", "T", "C", "G"))) return ("dna")
  
  if(all(unique_chars %in% c("A", "U", "C", "G"))) return ("rna")
  
  return ("aa")
}


fix_dataset_names <- function(Dataset, alias_csv){
  aliases <- read.csv(alias_csv, header = TRUE)
  for (i in seq_along(aliases$original)){
    Dataset <- ifelse(Dataset == aliases$original[i], aliases$replacement[i], Dataset)
  }
  return(Dataset)
}

fix_method_names <- function(Method, alias_csv){
  aliases <- read.csv(alias_csv, header = TRUE)
  for (i in seq_along(aliases$original)){
    Method <- ifelse(Method == aliases$original[i], aliases$replacement[i], Method)
  }
  return(factor(Method, 
                levels = aliases$replacement[aliases$order],
                labels = aliases$replacement[aliases$order],
                ordered = TRUE))
}

calc_stats2 <- function(mydf, refmethod = "EpitopeTransfer", blockEmbedder = FALSE, posthoc = "Dunnett"){
  
  mydf <- mydf %>% 
    group_by(Method) %>%
    arrange(Metric, Dataset) %>%
    ungroup()
  
  metrics <- unique(mydf$Metric)
  
  res <- vector("list", length(metrics))
  names(res) <- metrics
  for (i in seq_along(metrics)){
    tmp <- mydf %>% filter(Metric == metrics[i])
    
    if (blockEmbedder){
      mm <- aov(Value ~ Method + Embedder + Dataset, data = tmp)
    } else {
      mm <- aov(Value ~ Method + Dataset, data = tmp)
    }
    
    mymht <- glht(mm, linfct = mcp("Method" = posthoc))
    
    res[[i]]$df        <- tmp
    res[[i]]$aov.model <- mm
    res[[i]]$aov.summ  <- aov.summ <- summary(mm)
    res[[i]]$mht.summ  <- mht.summ <- summary(mymht)
    res[[i]]$mht.ci    <- mht.ci   <- confint(mymht)
    res[[i]]$summary   <- data.frame(Reference  = rep(refmethod, nrow(mht.ci$confint)),
                                    Challenger = gsub(paste0("- ", refmethod), "", unname(dimnames(mht.summ$linfct)[[1]])),
                                    Metric     = unique(mydf$Metric)[i],
                                    pValue     = mht.summ$test$pvalues,
                                    se         = mht.summ$test$sigma,
                                    diff       = - mht.ci$confint[, 1],
                                    ci.lower   = - mht.ci$confint[, 3],
                                    ci.upper   = - mht.ci$confint[, 2],
                                    p.str      = paste0("p = ", format(signif(mht.summ$test$pvalues, 3),
                                                                       scientific = TRUE)),
                                    ci.str     = sprintf("Mean of paired diffs: %2.3f (%2.3f, %2.3f)", 
                                                         - mht.ci$confint[, 1], 
                                                         - mht.ci$confint[, 3], 
                                                         - mht.ci$confint[, 2]))
  }
  
  return(res)
}






# calc_stats <- function(mydf, refmethod = "EpitopeTransfer"){
#   
#   mydf <- mydf %>% 
#     group_by(Method) %>%
#     arrange(Metric, Method, Dataset) %>%
#     ungroup()
#   
#   methods <- as.character(unique(mydf$Method))
#   methods <- methods[-which(methods == refmethod)]
#   lapply(seq_along(unique(mydf$Metric)),
#          function(i){
#            tmp <- mydf %>% filter(Metric == unique(mydf$Metric)[i])
#            x   <- tmp$Value[tmp$Method == refmethod]
#            res <- data.frame(Reference  = rep(refmethod, length(methods)),
#                              Challenger = NA_character_,
#                              Metric     = unique(mydf$Metric)[i],
#                              pValue = NA_real_,
#                              diff     = NA_real_,
#                              ci.lower = NA_real_,
#                              ci.upper = NA_real_)
#            for (j in seq_along(methods)){
#              y <- tmp$Value[tmp$Method == methods[j]]
#              cmp <- wilcox.test(x, y, paired = TRUE, conf.int = TRUE)
#              res$Challenger[j] <- methods[j]
#              res$pValue[j]     <- cmp$p.value
#              res$diff[j]       <- cmp$estimate
#              res$ci.lower[j]   <- cmp$conf.int[1]
#              res$ci.upper[j]   <- cmp$conf.int[2]
#            }
#            res$p.adj <- p.adjust(res$pValue, method = "fdr")
#            return(res)
#          }) %>%
#     bind_rows() %>%
#     mutate(p.adj.sci = paste0("p = ", format(signif(p.adj, 3),
#                                              scientific = TRUE)),
#            ci = sprintf("Median of paired diffs: %2.3f (%2.3f, %2.3f)", 
#                         diff, ci.lower, ci.upper))
# }
# 
# est_CIs <- function(mydf){
#   
#   mydf <- mydf %>% 
#     group_by(Method) %>%
#     arrange(Metric, Method, Dataset) %>%
#     ungroup()
#   
#   methods <- as.character(unique(mydf$Method))
#   lapply(seq_along(unique(mydf$Metric)),
#          function(i){
#            tmp <- mydf %>% filter(Metric == unique(mydf$Metric)[i])
#            res <- data.frame(Method   = rep(NA_real_, length(methods)),
#                              Metric   = unique(mydf$Metric)[i],
#                              median   = NA_real_,
#                              ci.lower = NA_real_,
#                              ci.upper = NA_real_)
#            for (j in seq_along(methods)){
#              x <- tmp$Value[tmp$Method == methods[j]]
#              cmp <- wilcox.test(x, conf.int = TRUE)
#              res$Method[j]   <- methods[j]
#              res$median[j]   <- cmp$estimate
#              res$ci.lower[j] <- cmp$conf.int[1]
#              res$ci.upper[j] <- cmp$conf.int[2]
#            }
#            return(res)
#          }) %>%
#     bind_rows() %>%
#     mutate(ci = sprintf("Median: %2.3f (%2.3f, %2.3f)", 
#                         median, ci.lower, ci.upper))
# }
# 
# my_auc10 <- function(truth, prob, posValue = "TRUE", negValue = "FALSE"){
#   # Implementation validated against pROC::roc
#   
#   if(sum(truth == posValue) == 0 || sum(truth == negValue) == 0) return(NA)
#   
#   tr  <- sort(unique(prob), decreasing = TRUE)
#   tpr <- fpr <- numeric(length(tr))
#   auc10 <- 0
#   for (i in seq_along(tr)){
#     pred <- ifelse(prob >= tr[i],
#                    yes = posValue,
#                    no  = negValue)
#     tpr[i] <- sum(pred == posValue & truth == posValue) / sum(truth == posValue)
#     fpr[i] <- sum(pred == posValue & truth == negValue) / sum(truth == negValue)
#     if (i > 1){
#       if(fpr[i] <= 0.1){
#         auc10 <- auc10 + (fpr[i] - fpr[i - 1]) * (tpr[i] + tpr[i - 1]) / 2
#       } else {
#         k = (.1 - fpr[i-1]) / (fpr[i] - fpr[i-1] + 1e-12)
#         tpr[i] <- tpr[i-1] + k * (tpr[i] - tpr[i-1])
#         fpr[i] <- 0.1
#         auc10 <- auc10 + (fpr[i] - fpr[i - 1]) * (tpr[i] + tpr[i - 1]) / 2
#         return(auc10)
#       }
#     }
#   }
#   return(auc10)
# }



# 
# my_f1 <- function(truth, pred, posValue = "TRUE", negValue = "FALSE"){
#   
#   TN <- as.numeric(sum(truth == negValue & pred == negValue))
#   FN <- as.numeric(sum(truth == posValue & pred == negValue))
#   FP <- as.numeric(sum(truth == negValue & pred == posValue))
#   TP <- as.numeric(sum(truth == posValue & pred == posValue))
#   
#   return(2 * TP / (2 * TP + FP + FN + 1e-16))
#   
# }
# 
# my_mcc <- function(truth, pred, posValue = "TRUE", negValue = "FALSE"){
#   
#   TN <- as.numeric(sum(truth == negValue & pred == negValue))
#   FN <- as.numeric(sum(truth == posValue & pred == negValue))
#   FP <- as.numeric(sum(truth == negValue & pred == posValue))
#   TP <- as.numeric(sum(truth == posValue & pred == posValue))
#   
#   mccNum  <- TP * TN - FP * FN
#   mccDen  <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
#   
#   return(mccNum / (mccDen + 1e-16))
#   
# }
# 
# my_auc <- function(truth, prob, posValue = "TRUE", negValue = "FALSE"){
#   
# myroc <- pROC::roc(truth, prob, levels = c(negValue, posValue))
# 
# return(as.numeric(pROC::auc(myroc)))
#   
# }


