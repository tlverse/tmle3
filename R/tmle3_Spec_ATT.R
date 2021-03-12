#' Defines a TML Estimator (except for the data)
#'
#'
#' @importFrom R6 R6Class
#'
#' @export
#
tmle3_Spec_ATT <- R6Class(
  classname = "tmle3_Spec_ATT",
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
    make_params = function(tmle_task, likelihood) {
      treatment_value <- self$options$treatment_level
      control_value <- self$options$control_level
      A_levels <- tmle_task$npsem[["A"]]$variable_type$levels
      if (!is.null(A_levels)) {
        treatment_value <- factor(treatment_value, levels = A_levels)
        control_value <- factor(control_value, levels = A_levels)
      }
      treatment <- define_lf(LF_static, "A", value = treatment_value)
      control <- define_lf(LF_static, "A", value = control_value)
      ate <- Param_ATT$new(likelihood, treatment, control)
      tmle_params <- list(ate)
      return(tmle_params)
    }
  ),
  active = list(),
  private = list()
)

#' All Treatment Specific Means
#'
#' O=(W,A,Y)
#' W=Covariates
#' A=Treatment (binary or categorical)
#' Y=Outcome (binary or bounded continuous)
#' @importFrom sl3 make_learner Lrnr_mean
#' @param treatment_level the level of A that corresponds to treatment
#' @param control_level the level of A that corresponds to a control or reference level
#' @export
tmle_ATT <- function(treatment_level, control_level) {
  # TODO: unclear why this has to be in a factory function
  tmle3_Spec_ATT$new(treatment_level, control_level)
}
