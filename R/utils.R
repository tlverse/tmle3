#' Wald-Style Confidence Intervals
#'
#' @importFrom stats qnorm
#'
#' @keywords internal
#
wald_ci <- function(est, se, alpha = 0.95) {
  z <- abs(stats::qnorm(p = (1 - alpha) / 2))
  ci_low <- est - z * se
  ci_high <- est + z * se
  return(cbind(ci_low, ci_high))
}

#' Summarize Estimates
#'
#' Generates a \code{data.table} summarizing results with inference
#'
#' @importFrom stats var
#'
#' @param task \code{tmle3_Task} containing the observed data of interest; the
#'  same as that passed to ..
#' @param estimates \code{list}, TMLE estimates of parameter and ICs from
#'  \code{\link{tmle3_Fit}$estimates}
#' @param param_types the types of the parameters being estimated
#' @param param_names the names of the parameters being estimated
#' @param init_psi the names of the parameters being estimated
#'
#' @return \code{data.table} summarizing results
#'
#' @export
#
summary_from_estimates <- function(task, estimates, param_types = NULL,
                                   param_names = NULL, init_psi = NULL) {
  psi <- unlist(lapply(estimates, `[[`, "psi"))
  IC <- lapply(estimates, `[[`, "IC")
  # for repeated measures, average IC values to get subject-level IC values
  if (length(unique(task$id)) < length(task$id)) {
    IC <- lapply(IC, function(x) {
      as.matrix(by(as.numeric(unlist(x)), as.numeric(task$id), mean))
    })
  }
  var_D <- unlist(lapply(IC, var))
  n <- sapply(IC, length)
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
  summary_dt <- as.data.table(list(param_types,
    param_names, init_psi, psi, se, ci,
    psi_transformed, ci_transformed
  ))
  setnames(summary_dt, c(
    "type", "param", "init_est", "tmle_est", "se", "lower",
    "upper", "psi_transformed", "lower_transformed",
    "upper_transformed"
  ))
  return(summary_dt)
}
