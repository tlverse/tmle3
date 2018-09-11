#' Defines a TML Estimator for the Risk Ratio
#'
#' Current limitations:
#' pretty much tailored to Param_TSM
#' see TODOs for places generalization can be added
#'
#' @importFrom R6 R6Class
#'
#' @export
#
tmle3_Spec_RR <- R6Class(
  classname = "tmle3_Spec_RR",
  portable = TRUE,
  class = TRUE,
  inherit = tmle3_Spec,
  public = list(
    initialize = function(baseline_level = 0, contrast = 1, ...) {
      # TODO: use sl3 param grabbing code
      options <- list(
        baseline_level = baseline_level,
        contrast_level = contrast
      )
      do.call(super$initialize, options)
    },
    make_params = function(tmle_task, likelihood) {
      baseline_level <- self$options$baseline_level
      contrast_level <- self$options$contrast_level

      intervention_base <- define_lf(LF_static, "A", value = baseline_level)
      intervention_cont <- define_lf(LF_static, "A", value = contrast_level)

      tsm_base <- Param_TSM$new(likelihood, intervention_base)
      tsm_cont <- Param_TSM$new(likelihood, intervention_cont)
      rr <- Param_delta$new(
        likelihood, delta_param_RR,
        list(tsm_base, tsm_cont)
      )
      tmle_params <- list(tsm_base, tsm_cont, rr)

      return(tmle_params)
    }
  ),
  active = list(),
  private = list()
)

#' Risk Ratio
#'
#' O = (W, A, Y)
#' W = Covariates
#' A = Treatment (binary or categorical)
#' Y = Outcome (binary or bounded continuous)
#'
#' @importFrom sl3 make_learner Lrnr_mean
#'
#' @param baseline_level The baseline risk group.
#' @param contrast_level The contrast risk group.
#'
#' @export
#
tmle_RR <- function(baseline_level, contrast_level) {
  # TODO: unclear why this has to be in a factory function
  tmle3_Spec_RR$new(baseline_level, contrast_level)
}
