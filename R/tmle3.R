#' TMLE from a tmle3_Spec object
#'
#' Using a tmle3_Spec object, fit a TMLE
#'
#' @param tmle_spec \code{\link{tmle3_Spec}}, defines the TMLE
#' @param data \code{data.frame}, the raw data
#' @param node_list \code{list}, defines which variables are which nodes
#' @param learner_list \code{list}, defines which learners are used to fit which likelihood factors
#' @return A \code{\link{tmle3_Fit}} object
#' @import data.table
#' @export
tmle3 <- function(tmle_spec, data, node_list, learner_list = NULL) {
  start_time <- proc.time()

  tmle_task <- tmle_spec$make_tmle_task(data, node_list)
  task_time <- proc.time()

  initial_likelihood <- tmle_spec$make_initial_likelihood(tmle_task, learner_list)
  likelihood_time <- proc.time()

  updater <- tmle_spec$make_updater()
  targeted_likelihood <- tmle_spec$make_targeted_likelihood(initial_likelihood, updater)

  tmle_params <- tmle_spec$make_params(tmle_task, targeted_likelihood)
  updater$tmle_params <- tmle_params
  params_time <- proc.time()

  fit <- fit_tmle3(tmle_task, targeted_likelihood, tmle_params, updater)
  fit_time <- proc.time()

  fit$set_timings(start_time, task_time, likelihood_time, params_time, fit_time)

  return(fit)
}
