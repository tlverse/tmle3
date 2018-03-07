#' TMLE from a tmle3_Spec object
#'
#' Using a tmle3_Spec object, fit a TMLE
#'
#' @param tmle_spec \code{\link{tmle3_Spec}}, defines the TMLE
#' @param data \code{data.frame}, the raw data
#' @param node_list \code{list}, defines which variables are which nodes
#' @param learner_list \code{list}, defines which learners are used to fit which likelihood factors
#' @return A \code{\link{tmle3_Fit}} object
#'
#' @export
tmle3 <- function(tmle_spec, data, node_list, learner_list = NULL) {
  start_time <- proc.time()

  tmle_task <- tmle_spec$make_tmle_task(data, node_list)
  task_time <- proc.time()

  likelihood <- tmle_spec$make_likelihood(tmle_task, learner_list)
  likelihood_time <- proc.time()

  tmle_params <- tmle_spec$make_params(tmle_task, likelihood)
  params_time <- proc.time()

  delta_params <- tmle_spec$make_delta_params()
  updater <- tmle_spec$make_updater(likelihood, tmle_params)
  fit <- fit_tmle3(tmle_task, likelihood, tmle_params, updater, delta_params)
  fit_time <- proc.time()

  fit$set_timings(start_time, task_time, likelihood_time, params_time, fit_time)

  return(fit)
}
