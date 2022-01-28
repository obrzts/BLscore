source("R/calculate_scores.R")
library(data.table)
library(dplyr)

clusterize_TCR <- function(sequence_df, id_col, chains, tmp_folder, scores_filename=NA, threshold=NA, ncores=1){
  
  ### INPUT CHECK

  # check chains argument
  if (!chains %in% c("AB", "B")){
    stop('Invalid chains argument (must be "AB" or "B")')
  }
  
  
  # check if all the needed columns are present
  beta_cols = c("v_beta", "j_beta", "cdr3_beta")
  missing_beta = beta_cols[!beta_cols %in% colnames(sequence_df)]
  if (length(missing_beta) > 0){
    stop(paste0("Coulums ", paste(missing_beta, collapse = ", "), " are missing in sequence_df"))
  }
  if (chains == "AB"){
    alpha_cols = c("v_alpha", "j_alpha", "cdr3_alpha")
    missing_alpha = alpha_cols[!alpha_cols %in% colnames(sequence_df)]
    if (length(missing_alpha) > 0){
      stop(paste0("Chains AB specified, but coulums ", 
                  paste(missing_alpha, collapse = ", "), 
                  " are missing in sequence_df"))
    }
  }
  
  # convert to data.table
  sequence_dt <- data.table::as.data.table(sequence_df)
  
  
  # filter out sequences with non-IMGT genes
  genes <- unique(c(sequence_dt$v_alpha, sequence_dt$v_beta, 
                    sequence_dt$j_alpha, sequence_dt$j_beta))
  unkn_genes =  genes[!genes %in% imgt_gene_list]
  if (length(unkn_genes) > 0){
    warning(paste("Unknown genes detected: ", paste(unkn_genes, collapse = ", "),
                  "Sequences carrying these genes will not be processed"))
    
    sequence_dt = sequence_dt[!sequence_dt$v_beta %in% unkn_genes]
    sequence_dt = sequence_dt[!sequence_dt$j_beta %in% unkn_genes]
    
    if (chains == "AB"){
      sequence_dt = sequence_dt[!sequence_dt$v_alpha %in% unkn_genes]
      sequence_dt = sequence_dt[!sequence_dt$j_alpha %in% unkn_genes]
    }
  }

  
  # filter out short CDR3 sequences
  n_short_beta = sum(nchar(sequence_dt$cdr3_beta) < 5)
  if (n_short_beta > 0){
    warning(paste0(n_short_beta, " sequences having short CDR3 beta (<5 aa) will not be processed"))
  }
  sequence_dt = sequence_dt[nchar(sequence_dt$cdr3_beta) > 4]
  
  if (chains == "AB"){
    n_short_alpha = sum(nchar(sequence_dt$cdr3_alpha) < 5)
    if (n_short_alpha > 0){
      warning(paste0(n_short_alpha, " sequences having short CDR3 alpha (<5 aa) will not be processed"))
    }
    sequence_dt = sequence_dt[nchar(sequence_dt$cdr3_alpha) > 4]
  }
  
  
  # filter out sequences with stop-codons
  has_stop_beta = grep("\\*", sequence_dt$cdr3_beta)
  if (chains == "AB"){
    has_stop_alpha = grep("\\*", sequence_dt$cdr3_alpha)
    has_stop = unique(has_stop_alpha, has_stop_beta)
  } else {
    has_stop = has_stop_beta
  }
  
  if (length(has_stop) > 0){
    warning(paste0(length(has_stop), " sequences having stop codon in CDR3 will not be processed"))
    sequence_dt = sequence_dt[-has_stop,]
  }
  
  
  # select or create unique receptor_id column and index the data.table
  if (missing(id_col)){
    sequence_dt[, receptor_id := 1:.N]
  
  } else {
    # check if id_col exists
    if (id_col %in% colnames(sequence_dt)){
      # check if id_col contains unique values
      if (nrow(sequence_dt) != length(unique(sequence_dt[[id_col]]))){
        warning(paste0("Provided id_col (", id_col, 
                       ") contains duplicated values, ids will be generated based on the order in the input table" ))
        sequence_dt[, receptor_id := 1:.N]
      } else {
        sequence_dt[, receptor_id := get(id_col)]
      }
    } else {
      warning(paste0("Column with name ", id_col, 
                     " is not found, ids will be generated based on the order in the input table" ))
      sequence_dt[, receptor_id := 1:.N]
    }
  }
  
  data.table::setkey(sequence_dt, receptor_id)
  
  
  ### CALCULATE SCORES
  scored_rec_pairs = get_scores(sequence_dt, chains, tmp_folder, scores_filename, ncores)
  
  
  ### DEFINE CLUSTERS
  if (is.na(threshold)){
    threshold = ifelse(chains == "AB", threshold_AB, threshold_B)
  }
  
  sim_rec_pairs = scored_rec_pairs[score > threshold, .(from_receptor_id, to_receptor_id, weight = score)]
  sim_rec_pairs$from_receptor_id <- as.character(sim_rec_pairs$from_receptor_id) # node ids should be characters otherwise messed up by igraph
  sim_rec_pairs$to_receptor_id <- as.character(sim_rec_pairs$to_receptor_id)
  
  g <- igraph::graph_from_data_frame(sim_rec_pairs,
                                     directed = F,
                                     vertices = as.character(sequence_df$receptor_id))
  clusters <- igraph::clusters(g)$membership
  
  clusters <- clusters %>%
    stack %>%
    dplyr::rename(cluster_id = values,
                  receptor_id = ind) %>%
    dplyr::mutate(receptor_id = as.integer(as.character(receptor_id))) %>%
    merge(sequence_dt, by = "receptor_id")
  
  return(clusters)
}


sequence_df = readRDS("test/test_B_big.Rds")
clusters <- clusterize_TCR(sequence_df, id_col="receptor_id", chains="B", tmp_folder=".", scores_filename="tmp.Rds", threshold=NA, ncores=7)

sequence_df = readRDS("test/test_AB_big.Rds")
clusters <- clusterize_TCR(sequence_df, id_col="receptor_id", chains="AB", tmp_folder=".", scores_filename="tmp.Rds", threshold=NA, ncores=7)

