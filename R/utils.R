#' @importFrom stats qnorm
wald_ci <- function(est, se, alpha=0.95) {
  z <- qnorm((1 + alpha) / 2)
  lower <- est - z * se
  upper <- est + z * se
  ci <- cbind(lower, upper)
  return(ci)
}

#' Summarize Estimates
#'
#' Generates a \code{data.table} summarizing results with inference
#' @importFrom stats var
#' @param estimates \code{list}, TMLE estimates of parameter and ICs from \code{\link{tmle3_Fit}$estimates}
#' @param param_types the types of the parameters we are estimating
#' @param param_names the names of the parameters we are estimating
#' @param init_psi the initial estimates
#' @return \code{data.table} summarizing results
#' @export
summary_from_estimates <- function(estimates, param_types = NULL, param_names = NULL, init_psi = NULL) {
  psi <- sapply(estimates, `[[`, "psi")
  IC <- sapply(estimates, `[[`, "IC")
  var_D <- apply(IC, 2, var)
  n <- length(estimates[[1]]$IC)
  se <- sqrt(var_D / n)
  ci <- wald_ci(psi, se)

  if (is.null(param_types)) {
    param_types <- rep(as.character(NA), length(estimates))
  }

  if (is.null(param_names)) {
    param_names <- rep(as.character(NA), length(estimates))
  }
  if (is.null(init_psi)) {
    init_psi <- rep(as.double(NA), length(estimates))
  }

  transforms <- sapply(estimates, function(estimate) {
    if (!is.null(estimate$transform)) {
      return(estimate$transform)
    } else {
      return(identity)
    }
  })
  apply_transform <- function(x, transform) {
    transform(x)
  }

  psi_transformed <- mapply(apply_transform, psi, transforms)
  ci_transformed <- mapply(apply_transform, ci, transforms)
  ci_transformed <- matrix(ci_transformed, nrow = nrow(ci), ncol = ncol(ci))
  summary_dt <- as.data.table(list(param_types, param_names, init_psi, psi, se, ci, psi_transformed, ci_transformed))
  setnames(summary_dt, c("type", "param", "init_est", "tmle_est", "se", "lower", "upper", "psi_transformed", "lower_transformed", "upper_transformed"))
  return(summary_dt)
}
