#' Defining Counterfactuals
#'
#' @param intervention_list ...
#' @param name ...
#'
#' @export
#
define_cf <- function(intervention_list) {
  if (inherits(intervention_list, "LF_base")) {
    intervention_list <- list(intervention_list)
  }

  intervention_nodes <- sapply(intervention_list, `[[`, "name")
  names(intervention_list) <- intervention_nodes

  return(intervention_list)
}
