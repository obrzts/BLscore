#' Create column with full junction sequence
#'
#' @param sequence_df A data.frame containing CDR3 sequences and J genes.
#' @param chain_to_modify "A", "B" or "AB"
#' @param anchor_tbl Table with anchor residues pre-filtered for a particular species
#'
#' @return Same data.frame with additional column junction_<alpha/beta>
#' @noRd
cdr3_to_junction <- function(sequence_df, chain_to_modify, anchor_tbl){
  chain_name = if (chain_to_modify == "B") "beta" else "alpha"
  j_name = paste("j", chain_name, sep = "_")
  junc_name = paste("junction", chain_name, sep = "_")
  cdr3_name = paste("cdr3", chain_name, sep = "_")

  # check unknown J genes
  known_j_beta = dplyr::filter(anchor_tbl, chain == chain_to_modify)$j
  genes <- unique(sequence_df[[j_name]])
  unkn_genes <- genes[!genes %in% known_j_beta]

  # merge tables and modify junction sequence
  res <- sequence_df |>
    merge(anchor_tbl |>
            dplyr::filter(chain == chain_to_modify) |>
            dplyr::select(j, anchor) |>
            dplyr::rename_with(~ c(j_name, "anchor"), c(j, anchor)),
          all.x = T) |>
    dplyr::mutate(!!junc_name := ifelse(is.na(anchor), NA,
                                        paste0("C", get(cdr3_name), anchor))) |>
    dplyr::select(-anchor)

  # notify about unknown genes
  if (length(unkn_genes) > 0) {
    message(paste("Unknown genes detected: ", paste(unkn_genes, collapse = ", "),
                  ". Sequences carrying these genes will get NA values"))
  }

  return(res)
}


#' Add CDR3 anchor residues
#'
#' @description Construct full TCR junction sequence by adding conservative
#' anchor residues (C and F/W) to CDR3 sequence.
#' @param sequence_df A data.frame containing CDR3 sequences and J genes. Columns
#' must be named cdr3_beta and j_beta and/or cdr3_alpha and j_alpha.
#' @param chains Which chains to process. "B" for beta chain only, "A" for alpha,
#'"AB" for both.
#' @param species From which species are the sequences. Currently only human and
#' mouse are supported.
#'
#' @return Same data.frame as input with additional columns junction_beta and/or
#' jucntion_alpha
#' @export
#'
#' @examples
#' # make example table without anchor residues
#' df <- example_TCR_df
#' df$cdr3_beta <- substr(df$junction_beta, 2, nchar(df$junction_beta) - 2)
#'
#' # generate column with beta chain junction sequence
#' df <- add_cdr3_anchors(df, "B", species="human")
add_cdr3_anchors <- function(sequence_df, chains, species="human"){

  ### INPUT CHECK
  # check chains argument
  if (!chains %in% c("AB", "B", "A")) {
    stop('Invalid chains argument (must be "A", "B", or "AB")')
  }

  # check if all the needed columns are present
  if (chains %in% c("B", "AB")){
    beta_cols <- c("j_beta", "cdr3_beta")
    missing_beta <- beta_cols[!beta_cols %in% colnames(sequence_df)]
    if (length(missing_beta) > 0) {
      stop(paste0("Chains ", chains, " are specified, but coulumn(s) ",
                  paste(missing_beta, collapse = ", "),
                  " are missing in sequence_df"))
    }
  }
  if (chains %in% c("A", "AB")){
    alpha_cols <- c("j_alpha", "cdr3_alpha")
    missing_alpha <- alpha_cols[!alpha_cols %in% colnames(sequence_df)]
    if (length(missing_alpha) > 0) {
      stop(paste0("Chains ", chains, " are specified, but coulumn(s) ",
                  paste(missing_alpha, collapse = ", "),
                  " are missing in sequence_df"))
    }
  }

  # check species
  if (!species %in% c("human", "mouse")) {
    stop('Invalid species (only human or mouse are supported)')
  } else {
    anchor_tbl <- j_anchor |>
      dplyr::filter(sp == species)
  }

  ### ADDING RESIDUES
  if (chains %in% c("B", "AB")){
    sequence_df <- sequence_df |>
      dplyr::mutate(j_beta_raw = j_beta,
                    j_beta = sub("\\*.+$", "", j_beta_raw))
    sequence_df <- cdr3_to_junction(sequence_df, chain_to_modify = "B", anchor_tbl)
    sequence_df$j_beta = sequence_df$j_beta_raw
    sequence_df$j_beta_raw <- NULL
  }
  if (chains %in% c("A", "AB")){
    sequence_df <- sequence_df |>
      dplyr::mutate(j_alpha_raw = j_alpha,
                    j_alpha = sub("\\*.+$", "", j_alpha_raw))
    sequence_df <- cdr3_to_junction(sequence_df, chain_to_modify = "A", anchor_tbl)
    sequence_df$j_alpha = sequence_df$j_alpha_raw
    sequence_df$j_alpha_raw <- NULL
  }

  return(sequence_df)
}
