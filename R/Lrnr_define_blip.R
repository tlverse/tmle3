#' Learner Class for Blip Regressions
#'
#' @param counterfactual_0 A \code{Counterfactual} object defining one possible
#'  intervention
#' @param counterfactual_1 A \code{Counterfactual} object defining a second
#'  possible intervention
#' @param V A character vector listing the nodes on which to regress the blip
#'  binary blip for now
#'
#' @export
#
make_blip_chain_function <- function(counterfactual_0, counterfactual_1, V) {
  stop("not done yet")
  ate <- Param_ATE$new(counterfactual_0, counterfactual_1)
  # A-IPW transform is just ate_ic+psi
  blip_chain <- function(likelihood, task) {
    ests <- ate$estimates(likelihood, task)
    blip <- ests$IC + ests$psi
    covariates <- task$get_data(, V)
  }
  return(blip_chain)
}
