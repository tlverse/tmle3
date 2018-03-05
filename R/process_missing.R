#' Preprocess Data to Handle Missing Variables
#'
#' Process data to account for missingness in preparation for TMLE
#'
#' Because IPCW-TMLE is not yet implemented, currently all rows with
#' missing outcomes (Y) or treatments (A) will be dropped from the analysis.
#'
#' Then covariates will be processed as follows:
#' \enumerate{
#'  \item any covariate with more than \code{max_p_missing} missingness will be dropped
#'  \item indicators of missingness will be generated
#'  \item missing values will be median-imputed
#' }
#'
#' @param data, \code{data.table}, containing the missing variables
#' @param node_list, \code{list}, what nodes serve what purpose
#' @param max_p_missing, \code{numeric}, what proportion of missing is tolerable? Beyond that, the variable will be dropped from the analysis
#' @return \code{list} containing the following elements:
#'  \itemize{
#'    \item \code{data}, the updated dataset
#'    \item \code{node_list}, the updated list of nodes
#'    \item \code{n_dropped}, the number of observations dropped
#'    \item \code{dropped_cols}, the variables dropped due to excessive missingness
#'
#'  }
#' @importFrom stats median
#' @export
process_missing <- function(data, node_list, max_p_missing = 0.5) {
  # todo - do IPCW instead of dropping these
  drop_vars <- c(node_list$A, node_list$Y)
  drop_rows <- data[, apply(is.na(.SD), 1, any), .SDcols = drop_vars]
  filtered <- data[!drop_rows]
  n_dropped <- sum(drop_rows)

  impute_vars <- node_list$W
  p_missing <- sapply(filtered[, impute_vars, with = FALSE], function(x) mean(is.na(x)))
  to_drop <- names(p_missing[(max_p_missing < p_missing)])
  to_impute <- names(p_missing[(0 < p_missing) & (p_missing < max_p_missing)])
  no_missing <- names(p_missing[p_missing == 0])
  no_missing <- c(no_missing, drop_vars)
  missing_indicators <- filtered[, lapply(.SD, is.na), .SDcols = to_impute]
  missing_names <- sprintf("delta_%s", to_impute)
  setnames(missing_indicators, missing_names)
  impute_median <- function(x) {
    value <- median(as.numeric(x[!is.na(x)]))
    x[is.na(x)] <- value
    x
  }
  imputed <- filtered[, lapply(.SD, impute_median), .SDcols = to_impute]
  all_observed <- filtered[, no_missing, with = FALSE]

  processed <- cbind(all_observed, imputed, missing_indicators)
  updated_nodes <- lapply(node_list, function(node) {
    node_no_missing <- intersect(no_missing, node)
    node_imputed <- c(to_impute, missing_names)[to_impute %in% node]
    return(c(node_no_missing, node_imputed))
  })

  return(list(data = processed, node_list = updated_nodes, n_dropped = n_dropped, dropped_cols = to_drop))
}
