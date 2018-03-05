#' Get Inference for Differentiable Functions of Parameters
#'
#' Using the functional delta method, get estimates and inference for a smooth function of parameters
#'
#' @param estimates \code{list}, TMLE estimates of parameter and ICs from \code{\link{tmle3_Fit}$estimates}
#' @param f \code{function}, the function of parameter values
#' @param df \code{function}, the derivative of the function wrt the parameter values
#' @return \code{list}, delta method estimates of parameter and ICS
#' @export
delta_method <- function(estimates, f, df) {
  psis <- lapply(estimates, `[[`, "psi")
  ICs <- lapply(estimates, `[[`, "IC")
  psi <- f(psis)
  IC <- df(ICs)
  list(psi = psi, IC = IC)
}

#' Linear Contrast Function
#'
#'  \eqn{X_2-X_1}
#'
#' For use with \code{\link{delta_method}}
#'
#' @param x numeric vector, containing two parameter estimates
#' @return \code{x[[2]] - x[[1]]}
#' @export
f_contrast <- function(x) {
  x[[2]] - x[[1]]
}


# todo: integrate with tmle_fit methods
#' Summarize Estimate
#'
#' Generates a \code{data.table} summarizing results with inference
#'
#' @param estimate \code{list}, TMLE estimates of parameter and ICs from \code{\link{tmle3_Fit}$estimates}
#' @return \code{data.table} summarizing results
#' @export
summary_from_estimate <- function(estimate) {
  ED2 <- mean(estimate$IC ^ 2)
  se <- sqrt(ED2) / sqrt(length(estimate$IC))
  ci <- wald_ci(estimate$psi, se)
  summary_dt <- data.table(estimate$psi, se, ci)
  setnames(summary_dt, c("tmle_est", "se", "lower", "upper"))

  return(summary_dt)
}
