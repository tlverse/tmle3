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
