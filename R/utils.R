#' Wald-Style Confidence Intervals
#'
#' @importFrom stats qnorm
#'
#' @keywords internal
#
wald_ci <- function(est, se, level = 0.95, q = NULL) {
  if (is.null(q)) {
    q <- abs(stats::qnorm(p = (1 - level) / 2))
  }

  ci_low <- est - q * se
  ci_high <- est + q * se
  return(cbind(ci_low, ci_high))
}

#' Summarize Estimates
#'
#' Generates a \code{data.table} summarizing results with inference
#'
#' @importFrom stats var cov
#' @importFrom mvtnorm qmvnorm
#'
#' @param task \code{tmle3_Task} containing the observed data of interest; the
#'  same as that passed to ..
#' @param estimates \code{list}, TMLE estimates of parameter and ICs from
#'  \code{\link{tmle3_Fit}$estimates}
#' @param param_types the types of the parameters being estimated
#' @param param_names the names of the parameters being estimated
#' @param init_psi the names of the parameters being estimated
#' @param simultaneous_ci if TRUE, calculate simulatenous confidence intervals
#'
#' @return \code{data.table} summarizing results
#'
#' @export
#
summary_from_estimates <- function(task, estimates, param_types = NULL,
                                   param_names = NULL, init_psi = NULL,
                                   simultaneous_ci = FALSE) {
  psi <- unlist(lapply(estimates, `[[`, "psi"))

  IC <- lapply(estimates, `[[`, "IC")
  IC <- do.call(cbind, IC)
  # for repeated measures, average IC values to get subject-level IC values
  if (length(unique(task$id)) < length(task$id)) {
    combined <- (by(IC, as.numeric(task$id), colMeans, simplify = FALSE))
    IC <- do.call(rbind, combined)
  }

  var_D <- cov(IC)
  n <- nrow(IC)
  se <- sqrt(diag(var_D) / n)
  level <- 0.95

  if (simultaneous_ci && (ncol(IC) > 1)) {
    rho_D <- var_D / sqrt(tcrossprod(diag(var_D)))
    q <- qmvnorm(level, tail = "both", corr = rho_D)$quantile
  } else {
    q <- abs(stats::qnorm(p = (1 - level) / 2))
  }

  ci <- wald_ci(psi, se, q = q)


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

  psi_lengths <- sapply(lapply(estimates, `[[`, "psi"), length)
  index_vec <- rep(seq_along(psi_lengths), psi_lengths)

  psi_transformed <- mapply(apply_transform, psi, transforms[index_vec])
  ci_transformed <- mapply(apply_transform, ci, transforms[index_vec])
  ci_transformed <- matrix(ci_transformed, nrow = nrow(ci), ncol = ncol(ci))
  summary_dt <- as.data.table(list(
    param_types[index_vec],
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

#' Get Empirical Mean of EIFs from Estimates
#' @param estimates a list of estimates objects
ED_from_estimates <- function(estimates) {
  IC <- lapply(estimates, `[[`, "IC")
  IC <- do.call(cbind, IC)
  ED <- colMeans(IC)

  return(ED)
}

#' Bound (Truncate) Likelihoods
#' @param x the likelihood values to bound
#' @param bounds Either a length two vector of c(lower,upper) or a lower bound, where the upper is then 1 - lower
#' @export
bound <- function(x, bounds) {
  lower <- bounds[[1]]
  if (length(bounds) > 1) {
    upper <- bounds[[2]]
  } else {
    upper <- 1 - lower
  }
  pmin(pmax(x, lower), upper)
}

#' Manually Train Likelihood Factor
#' The internal training process for likelihood factors is somewhat obtuse, so this function
#' does the steps to manually train one, which is helpful if you want to use a likelihood factor
#' independently of a likelihood object
#' @param lf the likelihood factor to train
#' @param tmle_task the task to use for training
#' @export
train_lf <- function(lf, tmle_task) {
  lf_fit <- lf$delayed_train(tmle_task)
  if (inherits(lf_fit, "Delayed")) {
    lf_fit <- lf_fit$compute()
  }

  lf$train(tmle_task, lf_fit)

  return(lf)
}
