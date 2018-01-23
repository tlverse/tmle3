#' nonparamteric likelihood (for baseline covariates only)
#'
#' @importFrom R6 R6Class
#'
#' @export
#
LF_np <- R6Class(
  classname = "Lf_np",
  portable = TRUE,
  class = TRUE,
  inherit = LF_base,
  public = list(
    initialize = function(name, ...) {
      private$.name <- name
    },
    get_prediction = function(task) {
      stop("nothing to predict")
    },
    get_likelihood = function(task, only_observed = FALSE) {
      weights <- task$weights
      return(weights/sum(weights))

    }
  ),
  active = list(
  ),
  private = list(
    .name = NULL
  )
)
