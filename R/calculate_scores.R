#' Pairwise sequence alignment using pre-defined score matrix and gap opening and extension
#' weights
#'
#' @param seq_vector1 Vector of patterns for alignment
#' @param seq_vector2 Vector of subjects for alignment
#'
#' @return Vector of alignment scores
#' @noRd
get_alignment_scores <- function(seq_vector1, seq_vector2){

  Biostrings::pairwiseAlignment(Biostrings::AAStringSet(seq_vector1),
                                Biostrings::AAStringSet(seq_vector2),
                                substitutionMatrix = submat,
                                gapOpening = gap_open,
                                gapExtension = gap_ext,
                                scoreOnly = T)

}


get_scores_paired <- function(dt,
                              receptor1, receptor2_vec,
                              output_file_prefix){

  # generate table with receptor pairs
  col_names <- c("receptor_id", "cdr3_beta", "v_beta", "j_beta",
                "cdr3_alpha", "v_alpha", "j_alpha")
  receptor_pairs <- dt[.(receptor2_vec), ..col_names]
  colnames(receptor_pairs) <- paste("to", col_names, sep = "_")
  suppressWarnings(receptor_pairs[,  paste("from", col_names, sep = "_") :=
                                    as.list(dt[.(receptor1), ..col_names])])
  receptor_pairs <- receptor_pairs %>%
    merge(trv_scores %>% dplyr::rename_all(function(x) paste0(x, "_beta")), all.x=T,
          by = c("from_v_beta", "to_v_beta")) %>%
    merge(trj_scores %>% dplyr::rename_all(function(x) paste0(x, "_beta")), all.x=T,
          by = c("from_j_beta", "to_j_beta")) %>%
    merge(trv_scores %>% dplyr::rename_all(function(x) paste0(x, "_alpha")), all.x=T,
          by = c("from_v_alpha", "to_v_alpha")) %>%
    merge(trj_scores %>% dplyr::rename_all(function(x) paste0(x, "_alpha")), all.x=T,
          by = c("from_j_alpha", "to_j_alpha"))

  # align sequences
  receptor_pairs$cdr3_alpha_score <- get_alignment_scores(receptor_pairs$from_cdr3_alpha,
                                                          receptor_pairs$to_cdr3_alpha)
  receptor_pairs$cdr3_beta_score <- get_alignment_scores(receptor_pairs$from_cdr3_beta,
                                                         receptor_pairs$to_cdr3_beta)

  # write output in a temporary file
  filename <- paste0(output_file_prefix, receptor1, ".Rds")
  saveRDS(receptor_pairs, filename)
}


get_scores_beta <- function(dt, receptor1, receptor2_vec,
                            output_file_prefix){

  # generate table with receptor pairs
  col_names <- c("receptor_id", "cdr3_beta", "v_beta", "j_beta")
  receptor_pairs <- dt[.(receptor2_vec), ..col_names]
  colnames(receptor_pairs) <- paste("to", col_names, sep = "_")
  suppressWarnings(receptor_pairs[,  paste("from", col_names, sep = "_") :=
                   as.list(dt[.(receptor1), ..col_names])])
  receptor_pairs <- receptor_pairs %>%
    merge(trv_scores %>% dplyr::rename_all(function(x) paste0(x, "_beta")), all.x=T,
          by = c("from_v_beta", "to_v_beta")) %>%
    merge(trj_scores %>% dplyr::rename_all(function(x) paste0(x, "_beta")), all.x=T,
          by = c("from_j_beta", "to_j_beta"))

  # align sequences
  receptor_pairs$cdr3_beta_score <- get_alignment_scores(receptor_pairs$from_cdr3_beta,
                                                         receptor_pairs$to_cdr3_beta)

  # write output in a temporary file
  filename <- paste0(output_file_prefix, receptor1, ".Rds")
  saveRDS(receptor_pairs, filename)
}


#' Calculate the logistic scores for given TCR pairs
#'
#' @param data A data.table with alignment scores for each TCR pair
#' @param chains Which chains are available ("AB" or "B")
#'
#' @return A vector of BL-scores
#'
#' @noRd
calc_BL_score <- function(data, chains){
  # get feature names
  model_coeff <- if (chains == "AB") model_coeff_AB else model_coeff_B
  feature_names <- names(model_coeff)
  feature_names <- feature_names[-grep("Intercept", feature_names)]

  # calculate scores
  coeff_mat <- as.matrix(model_coeff)
  data_mat <- as.matrix(data[, ..feature_names])

  1 / (1 + exp(-(data_mat %*% coeff_mat[feature_names,] + model_coeff["(Intercept)"])))
}


#' Calculate BL-scores (BLOSUM-logistic) for each pair of TCRs in a dataset
#'
#' @param sequence_dt A data.table with TCR sequence data
#' @inheritParams clusterize_TCR
#' @return A data.table with id of first TCR, second TCR and the score
#' @noRd
calculate_scores <- function(sequence_dt, chains, tmp_folder, scores_filename, ncores){
  res <- tryCatch({

    # create temporary folder for files with individual receptor alignment results
    # individual results are stored and then merged to avoid RAM overconsumption
    tmp_folder_full <- paste0(tmp_folder, "/BL_score_tmp/")
    system(paste0("mkdir ", tmp_folder_full))

    # legacy: rename junction to cdr3
    sequence_dt[, cdr3_beta := junction_beta]
    if (chains == "AB") {
      sequence_dt[, cdr3_alpha := junction_alpha]
    }

    # get alignment scores
    a <- Sys.time()
    if (chains == "AB") {
      x <- parallel::mclapply(seq_len(nrow(sequence_dt) - 1),
                              function(i) get_scores_paired(sequence_dt,
                                                            receptor1 = sequence_dt$receptor_id[i],
                                                            receptor2_vec = sequence_dt$receptor_id[(i+1):nrow(sequence_dt)],
                                                            tmp_folder_full),
                              mc.cores = ncores)
    } else if (chains == "B") {
      x <- parallel::mclapply(seq_len(nrow(sequence_dt) - 1),
                              function(i) get_scores_beta(sequence_dt,
                                                          receptor1 = sequence_dt$receptor_id[i],
                                                          receptor2_vec = sequence_dt$receptor_id[(i+1):nrow(sequence_dt)],
                                                          tmp_folder_full),
                              mc.cores = ncores)
    } else {
      stop("Unrecognized chains argument (specify AB or B)")
    }
    a1 <- Sys.time()

    # create one merged file and delete temporary files
    merged_file <- paste0(tmp_folder, "/BL_scores.csv")

    # get column names from one of the files
    out_colnames <- colnames(readRDS(paste0(tmp_folder_full, sequence_dt$receptor_id[1], ".Rds")))
    write(paste(out_colnames, collapse = "\t"), merged_file)

    files <- list.files(tmp_folder_full, full.names = T)
    graph <- lapply(files, readRDS) %>%
      data.table::rbindlist()
    for (file in list.files(tmp_folder_full, full.names = T)) {
      system(paste0("rm -r ", file))
    }

    # load the file and calculate BL_score
    graph[, score := calc_BL_score(graph, chains)]
    cols <- c("from_receptor_id", "to_receptor_id", "score")
    graph <- graph[, ..cols]

    # save scores file if requested
    if (!is.na(scores_filename)) {
      if (grepl("[Rr]ds$", scores_filename)) {
        saveRDS(graph, scores_filename)
      } else {
        utils::write.table(graph, scores_filename, quote = F, sep = "\t", row.names = F)
      }
    }

    b <- Sys.time()
    time_diff <- b-a
    print(paste0("Scores were calculated in ", round(time_diff, 3), " ", units(time_diff)))

    time_diff <- b-a1
    print(paste0("Aligned in ", round(time_diff, 3), " ", units(time_diff)))

    graph
  },

  error = function(cond){
    message(cond)
    return(NULL)
  },

  warning = function(cond){
    message(cond)
    return(NULL)
  },

  finally = {
    # legacy: remove cdr3 column
    sequence_dt[, cdr3_beta := NULL]
    if (chains == "AB") {
      sequence_dt[, cdr3_alpha := NULL]
    }


    # remove temporary files folder
    tmp_folder_full = paste0(tmp_folder, "/BL_score_tmp")
    if (file.exists(tmp_folder_full)) {
      system(paste0("rm -r ", tmp_folder_full))
    }

    # remove temporary file with scores if exists
    if (file.exists(paste0(tmp_folder, "/BL_scores.csv"))) {
      system(paste0("rm -r ", tmp_folder, "/BL_scores.csv"))
    }
  })
  return(res)
}

