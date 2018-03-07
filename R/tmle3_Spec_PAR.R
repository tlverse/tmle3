#' Defines a tmle (minus the data)
#'
#' Current limitations:
#' pretty much tailored to Param_TSM
#' see todos for places generalization can be added
#' @importFrom R6 R6Class
#'
#' @export
#
tmle3_Spec_PAR <- R6Class(
  classname = "tmle3_Spec_PAR",
  portable = TRUE,
  class = TRUE,
  inherit=tmle3_Spec,
  public = list(
    initialize = function(baseline_level=1, ...) {
      #todo: use sl3 param grabbing code
      params <- list(baseline_level=baseline_level)
      do.call(super$initialize,params)
    },
    make_params = function(tmle_task, likelihood) {
      baseline_level <- self$params$baseline_level
      intervention <- define_lf(LF_static, "A", value = baseline_level)
      tsm <- Param_TSM$new(likelihood, intervention)
      mean_param <- Param_mean$new(likelihood)
      tmle_params <- list(mean_param, tsm)
      return(tmle_params)
    }
  ),
  active = list(
  ),
  private = list(
  )
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