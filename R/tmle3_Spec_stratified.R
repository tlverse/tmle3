#' Defines a Stratified TML Estimator (except for the data)
#'
#' @importFrom R6 R6Class
#' @importFrom assertthat assert_that
#'
#' @export
#
tmle3_Spec_stratified <- R6Class(
  classname = "tmle3_Spec_stratified",
  portable = TRUE,
  class = TRUE,
  inherit = tmle3_Spec,
  public = list(
    initialize = function(base_spec, base_estimate = TRUE, ...) {
      private$.base_spec <- base_spec
      private$.base_estimate <- base_estimate
    },
    make_tmle_task = function(data, node_list, ...) {
      private$.strata_variable = node_list$V
      # Initial estimate should include V as covariate
      # Note: Upon current structure, the easiest way is to add V into W.
      #       This won't cause problem in calculation, 
      #       but nodes dependency graph will be off.
      node_list$W = c(node_list$W, node_list$V)
      
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
      
      tmle_params <- strat_params
      if (private$.base_estimate) {
        tmle_params <- c(base_params, tmle_params)
      }
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
    .strata_variable = NULL,
    .base_estimate = NULL
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
#' @param base_estimate Indicate whether to report base parameter.
#'
#' @export
tmle_stratified <- function(base_spec, base_estimate = TRUE) {
  # TODO: unclear why this has to be in a factory function
  tmle3_Spec_stratified$new(base_spec, base_estimate = base_estimate)
}
