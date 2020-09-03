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
submodel_logit <- function(eps, offset, X) {


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
submodel_density <- function(eps, offset, X) {
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
  estimate <- bound(estimate, c(0.0005))
  preds <- -1 * log(estimate)
  return(preds)
}

#' Function that returns the necessary updating information/functions for
#' a given optimization type in the targetting step.
#' Currently supports targetting through logistic submodels or (1 + eps D) P likelihood submodels
#'
#' @param optim_type A targeting method for the optimization.
#'
#'
#' @export
#
submodel_spec = function(optim_type = c("logistic", "EIC")){
  type <- match.arg(optim_type)
  plug_f <- function(x, ...){ stop("GLM based optimization does not work for this submodel.")}
  if(type == "logistic"){
    spec <- list(submodel = submodel_logit, loss_function = logistic_loss, family = binomial(), offset_tranform = stats::qlogis)
  } else if(type == "EIC") {
    spec <- list(submodel = submodel_density, loss_function = loglik_loss, family = plug_f, offset_tranform = plug_f)
  }
  return(spec)
}

