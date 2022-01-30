#' TCR sequence clustering
#'
#' @description Find clusters of similar TCRs that are likely to recognize the
#' same epitope.
#' @param sequence_df A data.frame containing TCR sequence data. Each row must
#' describe a unique TCR sequence. Following fields are required:
#' * junction_beta - amino acid sequence of CDR3 plus the two flanking conserved residues
#' * v_beta, j_beta - V and J gene with or without allele; allele information is
#' not used for score calculation.
#'
#' If chains="AB" junction_alpha, v_alpha and j_alpha must be provided too.
#' @param id_col Name of a column with unique ids for each TCR (optional).
#' @param chains Which chains to cluster. "B" for beta chain only, "AB" for paired
#' alpha and beta chains.
#' @param tmp_folder Path to a directory for storing temporary files which are
#' deleted when clustering is finished.
#' @param scores_filename If a character string for naming a file is provided
#' BL-scores of each TCR pair will be exported to this file. Supported formats: .Rds, .csv.
#' @param threshold Clustering threshold (optional).
#' @param ncores The number of cores to use for parallel computation (default = 1).
#'
#' @return A data.frame containing same information as sequence_df plus the cluster ids.
#' If scores_filename is provided a file with pairwise BL-scores is created.
#' @export
#'
#' @details
#' The default clustering thresholds were defined to optimally detect clusters of
#' TCRs recognizing the same epitope. If instead of full junction only CDR3 sequence
#' witout flanking residues is provided the scores will be overestimated which may
#' lead to wrong cluster assignment.
#'
#' @examples
#' head(example_TCR_df)
#' clusterize_TCR(example_TCR_df, chains="AB", tmp_folder=".", ncores=4)
#'
clusterize_TCR <- function(sequence_df, chains, tmp_folder, id_col,
                           scores_filename=NA, threshold=NA, ncores=1){

  ### INPUT CHECK

  # check chains argument
  if (!chains %in% c("AB", "B")) {
    stop('Invalid chains argument (must be "AB" or "B")')
  }


  # check if all the needed columns are present
  beta_cols <- c("v_beta", "j_beta", "junction_beta")
  missing_beta <- beta_cols[!beta_cols %in% colnames(sequence_df)]
  if (length(missing_beta) > 0) {
    stop(paste0("Coulumn(s) ", paste(missing_beta, collapse = ", "),
                " are missing in sequence_df"))
  }
  if (chains == "AB") {
    alpha_cols <- c("v_alpha", "j_alpha", "junction_alpha")
    missing_alpha <- alpha_cols[!alpha_cols %in% colnames(sequence_df)]
    if (length(missing_alpha) > 0) {
      stop(paste0("Chains AB specified, but coulumn(s) ",
                  paste(missing_alpha, collapse = ", "),
                  " are missing in sequence_df"))
    }
  }

  # convert to data.table
  sequence_dt <- data.table::as.data.table(sequence_df)

  # remove alleles
  data.table::setnames(sequence_dt, old = c('v_beta','j_beta'), new = c('v_beta_raw','j_beta_raw'))
  sequence_dt[, v_beta := sub("\\*.+$", "", v_beta_raw)]
  sequence_dt[, j_beta := sub("\\*.+$", "", j_beta_raw)]
  if (chains == "AB") {
    data.table::setnames(sequence_dt, old = c('v_alpha','j_alpha'), new = c('v_alpha_raw','j_alpha_raw'))
    sequence_dt[, v_alpha := sub("\\*.+$", "", v_alpha_raw)]
    sequence_dt[, j_alpha := sub("\\*.+$", "", j_alpha_raw)]
  }

  # filter out sequences with non-IMGT genes
  genes <- unique(c(sequence_dt$v_alpha, sequence_dt$v_beta,
                    sequence_dt$j_alpha, sequence_dt$j_beta))
  unkn_genes <- genes[!genes %in% imgt_gene_list]
  if (length(unkn_genes) > 0) {
    warning(paste("Unknown genes detected: ", paste(unkn_genes, collapse = ", "),
                  "Sequences carrying these genes will not be processed"))

    sequence_dt <- sequence_dt[!sequence_dt$v_beta %in% unkn_genes]
    sequence_dt <- sequence_dt[!sequence_dt$j_beta %in% unkn_genes]

    if (chains == "AB") {
      sequence_dt <- sequence_dt[!sequence_dt$v_alpha %in% unkn_genes]
      sequence_dt <- sequence_dt[!sequence_dt$j_alpha %in% unkn_genes]
    }
  }


  # filter out short junction sequences
  n_short_beta <- sum(nchar(sequence_dt$junction_beta) < 5)
  if (n_short_beta > 0) {
    warning(paste0(n_short_beta, " sequences having short junction in beta chain (<5 aa) will not be processed"))
    sequence_dt <- sequence_dt[nchar(sequence_dt$junction_beta) > 4]
  }


  if (chains == "AB") {
    n_short_alpha <- sum(nchar(sequence_dt$junction_alpha) < 5)
    if (n_short_alpha > 0) {
      warning(paste0(n_short_alpha, " sequences having short junction in alpha chain (<5 aa) will not be processed"))
      sequence_dt <- sequence_dt[nchar(sequence_dt$junction_alpha) > 4]
    }
  }


  # filter out sequences with stop-codons
  has_stop_beta <- grep("\\*", sequence_dt$junction_beta)
  if (chains == "AB") {
    has_stop_alpha <- grep("\\*", sequence_dt$junction_alpha)
    has_stop <- unique(has_stop_alpha, has_stop_beta)
  } else {
    has_stop <- has_stop_beta
  }

  if (length(has_stop) > 0) {
    warning(paste0(length(has_stop), " sequences having stop codon in junction will not be processed"))
    sequence_dt <- sequence_dt[-has_stop,]
  }


  # select or create unique receptor_id column and index the data.table
  if (missing(id_col)) {
    sequence_dt[, receptor_id := seq_len(.N)]

  } else {
    # check if id_col exists
    if (id_col %in% colnames(sequence_dt)) {
      # check if id_col contains unique values
      if (nrow(sequence_dt) != length(unique(sequence_dt[[id_col]]))) {
        warning(paste0("Provided id_col (", id_col,
                       ") contains duplicated values, ids will be generated based on the order in the input table" ))
        sequence_dt[, receptor_id := seq_len(.N)]
      } else {
        sequence_dt[, receptor_id := get(id_col)]
      }
    } else {
      warning(paste0("Column with name ", id_col,
                     " is not found, ids will be generated based on the order in the input table" ))
      sequence_dt[, receptor_id := seq_len(.N)]
    }
  }

  # check if data table is not empty
  if (nrow(sequence_dt) < 2) {
    stop(paste0(nrow(sequence_dt), " sequences, nothing to cluster!"))
  }

  data.table::setkey(sequence_dt, receptor_id)


  ### CALCULATE SCORES
  scored_rec_pairs <- calculate_scores(sequence_dt, chains, tmp_folder, scores_filename, ncores)


  ### DEFINE CLUSTERS
  if (is.na(threshold)){
    threshold <- if (chains == "AB") threshold_AB else threshold_B
  }

  sim_rec_pairs <- scored_rec_pairs[score > threshold, .(from_receptor_id, to_receptor_id, weight = score)]
  # node ids should be characters otherwise messed up by igraph
  sim_rec_pairs$from_receptor_id <- as.character(sim_rec_pairs$from_receptor_id)
  sim_rec_pairs$to_receptor_id <- as.character(sim_rec_pairs$to_receptor_id)

  g <- igraph::graph_from_data_frame(sim_rec_pairs,
                                     directed = F,
                                     vertices = as.character(sequence_dt$receptor_id))
  clusters <- igraph::clusters(g)$membership

  # return columns with initial V gene names to the table
  sequence_dt[, c("v_beta", "j_beta") := NULL]
  data.table::setnames(sequence_dt, old = c('v_beta_raw','j_beta_raw'), new = c('v_beta','j_beta'))
  if (chains == "AB") {
    sequence_dt[, c("v_alpha", "j_alpha") := NULL]
    data.table::setnames(sequence_dt, old = c('v_alpha_raw','j_alpha_raw'), new = c('v_alpha','j_alpha'))
  }

  clusters <- clusters %>%
    stack %>%
    dplyr::rename(cluster_id = values,
                  receptor_id = ind) %>%
    dplyr::mutate(receptor_id = as.integer(as.character(receptor_id))) %>%
    # merge to the full data.table to add all the info
    merge(sequence_dt, by = "receptor_id")

  return(clusters)
}
