#' Defines a tmle (minus the data)
#'
#' Current limitations:
#' pretty much tailored to Param_TSM
#' see TODOs for places generalization can be added
#' @importFrom R6 R6Class
#'
#' @export
#
tmle3_Spec_PAR <- R6Class(
  classname = "tmle3_Spec_PAR",
  portable = TRUE,
  class = TRUE,
  inherit = tmle3_Spec,
  public = list(
    initialize = function(baseline_level = 1, ...) {
      # TODO: use sl3 param grabbing code
      options <- list(baseline_level = baseline_level)
      do.call(super$initialize, options)
    },
    make_params = function(tmle_task, likelihood) {
      baseline_level <- self$options$baseline_level
      intervention <- define_lf(LF_static, "A", value = baseline_level)
      tsm <- Param_TSM$new(likelihood, intervention)
      mean_param <- Param_mean$new(likelihood)
      par <- Param_delta$new(likelihood, delta_param_PAR, list(tsm, mean_param))
      paf <- Param_delta$new(likelihood, delta_param_PAF, list(tsm, mean_param))
      tmle_params <- list(tsm, mean_param, par, paf)

      return(tmle_params)
    }
  ),
  active = list(),
  private = list()
)

#' PAR and PAF
#'
#' O=(W,A,Y)
#' W=Covariates
#' A=Treatment (binary or categorical)
#' Y=Outcome (binary or bounded continuous)
#' @importFrom sl3 make_learner Lrnr_mean
#' @param baseline_level, the baseline risk group
#' @export
tmle_PAR <- function(baseline_level) {
  # todo: unclear why this has to be in a factory function
  tmle3_Spec_PAR$new(baseline_level)
}
