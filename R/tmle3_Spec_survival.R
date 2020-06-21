#' Defines a TML Estimator (except for the data)
#'
#'
#' @importFrom R6 R6Class
#'
#' @export
#

tmle3_Spec_survival <- R6Class(
  classname = "tmle3_Spec_survival",
  portable = TRUE,
  class = TRUE,
  inherit = tmle3_Spec,
  public = list(
    initialize = function(treatment_level, control_level, target_times = NULL, variable_types = NULL,...) {
      super$initialize(
        # TODO: check variable types
        # TODO: support multi-level treatments and etc
        treatment_level = treatment_level,
        control_level = control_level, 
        variable_types = variable_types, 
        target_times = target_times, 
        ...
      )
    },
    make_tmle_task = function(data, node_list, ...) {
      variable_types <- self$options$variable_types
      
      tmle_task <- survival_tx_task(data, node_list, survival_tx_npsem, variable_types)

      return(tmle_task)
    },

    make_initial_likelihood = function(tmle_task, learner_list = NULL) {
      likelihood <- survival_tx_likelihood(tmle_task, learner_list)
      return(likelihood)
    },

    make_params = function(tmle_task, likelihood) {
      treatment_value <- self$options$treatment_level
      control_value <- self$options$control_level

      treatment <- define_lf(LF_static, "A", value = treatment_value)
      control <- define_lf(LF_static, "A", value = control_value)

      # TODO: currently support treatment specific
      param_surv <- Param_survival$new(likelihood, treatment, 
                                       target_times = self$optiosn$target_times, 
                                       outcome_node = "N")
      tmle_params <- list(param_surv)
      return(tmle_params)
    }
  ),
  active = list(),
  private = list()
)

# TODO
#' @importFrom sl3 make_learner Lrnr_mean
#' @param treatment_level the level of A that corresponds to treatment
#' @param control_level the level of A that corresponds to a control or reference level
#' @export
# TODO: check variable types
tmle_survival <- function(treatment_level, control_level, target_times = NULL, variable_types = NULL) {
  tmle3_Spec_survival$new(treatment_level, control_level, target_times, variable_types)
}
