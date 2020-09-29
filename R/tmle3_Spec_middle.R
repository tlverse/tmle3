#' Defines a TML Estimator (except for the data)
#'
#' Longitudinal Mediation Targets
#'
#' @importFrom R6 R6Class
#'
#' @export
#
tmle3_Spec_middle <- R6Class(
  classname = "tmle3_Spec_middle",
  portable = TRUE,
  class = TRUE,
  inherit = tmle3_Spec,
  public = list(
    initialize = function(treatment_level, control_level, ...) {
      super$initialize(
        treatment_level = treatment_level,
        control_level = control_level, ...
      )
    },
    make_tmle_task = function(data, node_list, ...) {
      variable_types <- self$options$variable_types
      tmle_task <- middle_task(data, node_list, variable_types)
      return(tmle_task)
    },
    make_initial_likelihood = function(tmle_task, learner_list = NULL) {
      # produce trained likelihood when likelihood_def provided

      if (!is.null(self$options$likelihood_override)) {
        likelihood <- self$options$likelihood_override$train(tmle_task)
      } else {
        likelihood <- middle_likelihood(tmle_task, learner_list)  # see middle_helper
      }

      return(likelihood)
    },
    make_params = function(tmle_task, likelihood, if_projection = NULL) {
      temp_names <- names(tmle_task$npsem)
      loc_A <- grep("A", temp_names)
      # ZW todo: in future can be dynamic
      treatment_value <- self$options$treatment_level
      control_value <- self$options$control_level
      A_levels <- tmle_task$npsem[[ temp_names[loc_A[1]] ]]$variable_type$levels
      if (!is.null(A_levels)) {
        treatment_value <- factor(treatment_value, levels = A_levels)
        control_value <- factor(control_value, levels = A_levels)
      }
      # list of intervention nodes as LF_static objects
      treatment <- lapply(temp_names[loc_A], function(eachA) {
        define_lf(LF_static, eachA, value = treatment_value)
      })
      control <- lapply(temp_names[loc_A], function(eachA) {
        define_lf(LF_static, eachA, value = control_value)
      })
      names(treatment) <- names(control) <- temp_names[loc_A]
      # treatment <- define_lf(LF_static, "A", value = treatment_value)
      # control <- define_lf(LF_static, "A", value = control_value)
      if (is.null(if_projection)) {
        middle <- Param_middle$new(likelihood, treatment, control, outcome_node = last(temp_names))
      } else if (if_projection) {
        middle <- Param_middle_projection$new(likelihood, treatment, control, outcome_node = last(temp_names))
      }

      tmle_params <- list(middle)
      return(tmle_params)
    }
  ),
  active = list(),
  private = list()
)

#' Longitudinal Mediation Targets
#'
#' O=(L0, A1, R1, Z1, L1, (Y1), ..., An, Rn, Zn, Ln, Yn)
#' L0=Baseline covariates
#' A: Treatment (binary or categorical)
#' Z: Mediators
#' R: (time-varying) Covariates before mediator
#' L: (time-varying) Covariates after mediator
#' Y=Outcome (binary or bounded continuous)
#' @importFrom sl3 make_learner Lrnr_mean
#' @param treatment_level the level of A that corresponds to treatment
#' @param control_level the level of A that corresponds to a control or reference level
#' @export
tmle_middle <- function(treatment_level, control_level) {
  # TODO: unclear why this has to be in a factory function
  tmle3_Spec_middle$new(treatment_level, control_level)
}
