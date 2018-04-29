#' Defines a tmle (minus the data)
#'
#' Current limitations:
#' pretty much tailored to Param_TSM
#' see todos for places generalization can be added
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
    initialize = function(baseline_level=0, contrast=1,...) {
      # todo: use sl3 param grabbing code
      params <- list(baseline_level = baseline_level, contrast_level=contrast)
      do.call(super$initialize, params)
    },
    make_params = function(tmle_task, likelihood) {

      baseline_level <- self$params$baseline_level
      contrast_level <- self$params$contrast_level

      intervention_base <- define_lf(LF_static, "A", value = baseline_level)
      intervention_cont <- define_lf(LF_static, "A", value = contrast_level)

      tsm_base <- Param_TSM$new(likelihood, intervention_base)
      tsm_cont <- Param_TSM$new(likelihood, intervention_cont)

      tmle_params <- list(tsm_base,tsm_cont)
      return(tmle_params)
    },
    make_delta_params = function() {
      delta_params <- list(delta_param_RR)
    }
  ),
  active = list(),
  private = list()
)

#' Risk ratio
#'
#' O=(W,A,Y)
#' W=Covariates
#' A=Treatment (binary or categorical)
#' Y=Outcome (binary or bounded continuous)
#' @importFrom sl3 make_learner Lrnr_mean
#' @param baseline_level, the baseline risk group
#' @export
tmle_RR <- function(baseline_level, contrast_level) {
  # todo: unclear why this has to be in a factory function
  tmle3_Spec_RR$new(baseline_level, contrast_level)
}
