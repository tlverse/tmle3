#' @export
submodel_logit <- function(eps, X, offset){
  preds <- plogis(qlogis(offset) + X %*% eps)
  return(preds)
}
