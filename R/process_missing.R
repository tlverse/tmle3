impute_median <- function(x) {
  value <- median(as.numeric(x[!is.na(x)]))
  x[is.na(x)] <- value
  return(x)
}

#' @importFrom stats aggregate
impute_mode <- function(x) {
  count_df <- aggregate(count ~ x, data = data.frame(count = 1, x = x), sum)
  value <- count_df$x[which.max(count_df$count)]
  x[is.na(x)] <- value
  return(x)
}

impute_by_type <- function(x) {
  if (is.factor(x) || is.character(x)) {
    return(impute_mode(x))
  } else {
    return(impute_median(x))
  }
}


#' Preprocess Data to Handle Missing Variables
#'
#' Process data to account for missingness in preparation for TMLE
#'
#' Rows where there is missingness in any of the \code{complete_nodes} will be
#' dropped. Then, missingness will be median-imputed for the variables in the \code{impute_nodes}.
#' Indicator variables of missingness will be generated for these nodes.
#'
#' Then covariates will be processed as follows:
#' \enumerate{
#'  \item any covariate with more than \code{max_p_missing} missingness will be dropped
#'  \item indicators of missingness will be generated
#'  \item missing values will be median-imputed
#' }
#'
#' @param data, \code{data.table}, containing the missing variables
#' @param node_list, \code{list}, what variables comprise each node
#' @param complete_nodes, \code{character vector}, nodes we must observe
#' @param impute_nodes, \code{character vector}, nodes we will impute
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
process_missing <- function(data, node_list, complete_nodes = c("A", "Y"), impute_nodes = NULL, max_p_missing = 0.5) {
  data <- as.data.table(data)
  if (is.null(impute_nodes)) {
    impute_nodes <- setdiff(names(node_list), complete_nodes)
  }

  # drop rows where there is missingness for nodes that we are required to observe
  drop_vars <- unlist(node_list[complete_nodes])
  drop_rows <- data[, apply(is.na(.SD), 1, any), .SDcols = drop_vars]
  filtered <- data[!drop_rows]
  n_dropped <- sum(drop_rows)

  # median impute the other nodes and build indicators
  impute_vars <- unlist(node_list[impute_nodes])
  p_missing <- sapply(filtered[, impute_vars, with = FALSE], function(x) mean(is.na(x)))


  # nodes that are already complete
  no_missing <- names(p_missing[p_missing == 0])
  no_missing <- c(no_missing, drop_vars)
  processed <- filtered[, no_missing, with = FALSE]

  # nodes to impute
  to_impute <- names(p_missing[(0 < p_missing) & (p_missing < max_p_missing)])
  if (length(to_impute) > 0) {
    missing_indicators <- filtered[, lapply(.SD, function(x) as.numeric(!is.na(x))), .SDcols = to_impute]
    missing_names <- sprintf("delta_%s", to_impute)
    setnames(missing_indicators, missing_names)

    imputed <- filtered[, lapply(.SD, impute_by_type), .SDcols = to_impute]
    processed <- cbind(processed, imputed, missing_indicators)
  } else {
    missing_names <- c()
  }

  # nodes with too much missingness
  to_drop <- names(p_missing[(max_p_missing < p_missing)])

  updated_nodes <- lapply(node_list, function(node) {
    node_no_missing <- intersect(no_missing, node)
    node_imputed <- c(to_impute, missing_names)[to_impute %in% node]
    return(c(node_no_missing, node_imputed))
  })

  return(list(data = processed, node_list = updated_nodes, n_dropped = n_dropped, dropped_cols = to_drop))
}
