#' Defines a TML Estimator (except for the data)
#'
#'
#' @importFrom R6 R6Class
#'
#' @export
#
tmle3_Spec_Threshold <- R6Class(
  classname = "tmle3_Spec_Threshold",
  portable = TRUE,
  class = TRUE,
  inherit = tmle3_Spec,
  public = list(
    initialize = function(threshold_values = NULL, num_thresholds = 10, cdf_bins = 10, ...) {
      super$initialize(
        threshold_values = threshold_values,
        num_thresholds = num_thresholds,
        cdf_bins = cdf_bins, ...
      )
    },
    make_initial_likelihood = function(tmle_task, learner_list = NULL) {
      # produce trained likelihood when likelihood_def provided
      cdf_bins <- self$options$cdf_bins
      threshold_values <- self$options$threshold_values
      if(is.null(threshold_values)) {
        num_thresholds <- self$options$num_thresholds
        threshold_values <- unique(unlist(quantile(tmle_task$get_tmle_node("A", format = T)[[1]], seq(0.05, 0.95, length.out = num_thresholds ))))
        print(threshold_values)
        print(data.table(tmle_task$get_tmle_node("A", format = T)[[1]]))
        private$.options$threshold_values <- threshold_values
      }
      if (!is.null(self$options$likelihood_override)) {
        likelihood <- self$options$likelihood_override$train(tmle_task)
      } else {
        likelihood <- threshold_likelihood(tmle_task, learner_list, threshold_values, cdf_bins)
      }

      return(likelihood)
    },
    make_params = function(tmle_task, likelihood) {
      thres <- Param_thresh$new(likelihood, self$options$threshold_values)
      tmle_params <- list(thres)
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
tmle_Threshold <- function(threshold_values, num_thresholds, cdf_bins, ...) {
  # TODO: unclear why this has to be in a factory function
  tmle3_Spec_Threshold$new(threshold_values, num_thresholds, cdf_bins, ...)
}
