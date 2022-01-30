#' Sample TCR sequence data
#'
#' An exemplary dataset containing junction sequences and V and J genes of 1000
#' TCRs randomly sampled from VDJdb database.
#'
#' @format A data frame with 2000 rows and 6 variables:
#' \describe{
#'   \item{v_beta}{beta chain V gene}
#'   \item{j_beta}{beta chain J gene}
#'   \item{v_alpha}{alpha chain J gene}
#'   \item{j_alpha}{alpha chain J gene}
#'   \item{junction_alpha}{beta chain junction sequence}
#'   \item{junction_beta}{alpha chain junction sequence}
#'   \item{id}{unique id}
#' }
#' @source \url{https://vdjdb.cdr3.net/}
#' @references Dmitry V Bagaev, Renske M A Vroomans, Jerome Samir, Ulrik Stervbo,
#' Cristina Rius, Garry Dolton, Alexander Greenshields-Watson, Meriem Attaf,
#' Evgeny S Egorov, Ivan V Zvyagin, Nina Babel, David K Cole, Andrew J Godkin,
#' Andrew K Sewell, Can Kesmir, Dmitriy M Chudakov, Fabio Luciani, Mikhail Shugay,
#' VDJdb in 2019: database extension, new analysis infrastructure and a T-cell
#' receptor motif compendium. \emph{Nucleic Acids Research, gkz874}, 2019
"example_TCR_df"
