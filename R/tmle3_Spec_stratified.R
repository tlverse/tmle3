#' Defines a Stratified TML Estimator (except for the data)
#'
#' @importFrom R6 R6Class
#'
#' @export
#
tmle3_Spec_stratified <- R6Class(
  classname = "tmle3_Spec_stratified",
  portable = TRUE,
  class = TRUE,
  inherit = tmle3_Spec,
  public = list(
    initialize = function(base_spec, strata_variable, ...) {
      private$.base_spec <- base_spec
      private$.strata_variable <- strata_variable
    },
    make_tmle_task = function(data, node_list, ...) {
      tmle_task <- self$base_spec$make_tmle_task(data, node_list, ...)

      return(tmle_task)
    },
    make_initial_likelihood = function(tmle_task, learner_list = NULL) {
      initial_likelihood <- self$base_spec$make_initial_likelihood(
        tmle_task,
        learner_list
      )

      return(initial_likelihood)
    },
    make_updater = function(...) {
      updater <- self$base_spec$make_updater(...)
      return(updater)
    },
    make_targeted_likelihood = function(likelihood, updater) {
      targeted_likelihood <- self$base_spec$make_targeted_likelihood(
        likelihood,
        updater
      )

      return(targeted_likelihood)
    },
    make_params = function(tmle_task, targeted_likelihood) {
      base_params <- self$base_spec$make_params(tmle_task, targeted_likelihood)
      strat_params <- lapply(base_params, function(base_param) {
        define_param(
          Param_stratified, targeted_likelihood,
          base_param, self$strata_variable
        )
      })
      tmle_params <- c(base_params, strat_params)
      return(tmle_params)
    }
  ),
  active = list(
    base_spec = function() {
      return(private$.base_spec)
    },
    strata_variable = function() {
      return(private$.strata_variable)
    }
  ),
  private = list(
    .base_spec = NULL,
    .strata_variable = NULL
  )
)

#' Stratified version of TML estimator from other Spec classes
#'
#' O=(W,A,Y)
#' W=Covariates
#' A=Treatment (binary or categorical)
#' Y=Outcome (binary or bounded continuous)
#'
#' @importFrom sl3 make_learner Lrnr_mean
#'
#' @param base_spec An underlying spec to stratify.
#' @param strata_variable The variable(s) to use for stratification.
#'
#' @export
tmle_stratified <- function(base_spec, strata_variable) {
  # TODO: unclear why this has to be in a factory function
  tmle3_Spec_stratified$new(base_spec, strata_variable)
}
