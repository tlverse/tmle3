#' Logistic Submodel Fluctuation
#'
#' @param eps ...
#' @param X ...
#' @param offset ...
#'
#' @importFrom stats plogis qlogis
#'
#' @export
#
submodel_logit <- function(eps, X, offset) {
  preds <- stats::plogis(stats::qlogis(offset) + X %*% eps)
  return(preds)
}

#' Logistic Loss Function
#'
#' @param estimate ...
#' @param observed ...
#'
#'
#' @export
#
logistic_loss <- function(estimate, observed) {
  -1 * ifelse(observed == 1, log(estimate), log(1 - estimate))
}


#' Density Submodel Fluctuation
#' Standard (1+ eps D)P submodel through densities
#' @param eps ...
#' @param X ...
#' @param offset ...
#'
#' @importFrom stats plogis qlogis
#'
#' @export
#
submodel_density <- function(eps, X, offset) {
  preds <- (1 + X %*% eps) * offset
  return(preds)
}


#' Log-likelihood Loss Function
#'
#' @param estimate ...
#' @param observed ...
#'
#'
#' @export
#
loglik_loss <- function(estimate, observed) {
  preds <- -1 * log(estimate)
  return(preds)
}

