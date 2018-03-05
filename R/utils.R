#' @importFrom stats qnorm
wald_ci <- function(est, se, alpha=0.95) {
  z <- qnorm((1 + alpha) / 2)
  lower <- est - z * se
  upper <- est + z * se
  ci <- cbind(lower, upper)
  return(ci)
}
